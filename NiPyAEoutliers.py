#!/usr/bin/env python3
# NiPyAEoutliers.py
# Renxiang Qiu

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, ConcatDataset
from torch_geometric.loader import DataLoader as GeometricDataLoader
from torch_geometric.data import Data
from torch_geometric.utils import dense_to_sparse
import torch_geometric.nn as pyg_nn
import os
from scipy.io import savemat


# ===========================================
# 1) Retrieve inputs from MATLAB
# ===========================================
#  MATLAB calls:
#    learned_betas = pyrunfile("NiPyAEoutliers.py", "learned_betas",
#                              datain=beta_values, binatry_matrix=neighbourgh_matrix, saveGAE = saveGAE);
#

# ===========================================
datain = np.array(datain)          # shape: [nBeta, nTime, nChan, nSubject] (adjust as needed)
binatry_matrix = np.array(binatry_matrix)
saveGAE = str(saveGAE)            # "yes", "no", or "xai"

# Parse shapes
nBeta     = datain.shape[0]
nTime     = datain.shape[1]
nChan     = datain.shape[2]
nSubject  = datain.shape[3]

# Create folder to save models weights and XAI results
if saveGAE == "yes" or saveGAE == "xai":
    newdir = os.path.join(os.getcwd(), "Group_outlier_parametrization")
    os.makedirs(newdir, exist_ok=True)
    model_weights_dir = os.path.join(newdir, "model weights")
    os.makedirs(model_weights_dir, exist_ok=True)
if saveGAE == "xai":
    import Explainable_GAE as xGAE
    
    LP_dir = os.path.join(newdir, "linear probing")
    MK_dir = os.path.join(newdir, "masking")
    AB_dir = os.path.join(newdir, "ablation")
    AM_dir = os.path.join(newdir, "activation maximization")
    os.makedirs(LP_dir, exist_ok=True); os.makedirs(MK_dir, exist_ok=True)
    os.makedirs(AB_dir, exist_ok=True); os.makedirs(AM_dir, exist_ok=True)

# ===========================================
# 2) Set device & convert adjacency
# ===========================================
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
edge_index, _ = dense_to_sparse(torch.tensor(binatry_matrix, dtype=torch.float32))

# ===========================================
# 3) Dataset builder
# ===========================================
def build_concat_dataset_for_beta(datain, beta_idx):
    """
    Build a ConcatDataset across all subjects for one beta index.
    `datain` is [nBeta, nTime, nChan, nSubject].
    Returns: ConcatDataset of length nSubject.
    """
    subject_datasets = []
    for s in range(nSubject):
        # Extract slice for one subject and one beta
        # shape: [nTime, nChan]
        single_beta_data = datain[beta_idx, :, :, s]

        # If each node is a channel and its feature dimension is time,
        # we might want single_beta_data to be [nChan, nTime].
        # If so, transpose here:
        single_beta_data = single_beta_data.T  # shape: [nChan, nTime]

        # For an autoencoder, features=labels
        x_tensor = torch.tensor(single_beta_data, dtype=torch.float32)
        x_tensor = x_tensor.unsqueeze(0)

        subject_datasets.append(TensorDataset(x_tensor, x_tensor))

    return ConcatDataset(subject_datasets)

# ===========================================
# 4) PyG Graph Dataset
# ===========================================
class EEGGraphDataset(torch.utils.data.Dataset):
    """
    Wraps a standard dataset so that each item is a PyG Data object:
      x -> node features (channels)
      edge_index -> adjacency
    """
    def __init__(self, eeg_data, edge_index):
        self.eeg_data = eeg_data
        self.edge_index = edge_index

    def __len__(self):
        return len(self.eeg_data)

    def __getitem__(self, idx):
        # item: (x, x) because it's autoencoder
        sample = self.eeg_data[idx]
        x = sample[0]  # shape: [nChan, nTime] (if you used transpose above)
        # Convert x -> float32, build PyG Data
        graph_data = Data(x=torch.tensor(x, dtype=torch.float32),
                          edge_index=self.edge_index)
        return graph_data

# ===========================================
# 5) Define GraphAutoencoder
# ===========================================
class GraphAutoencoder(nn.Module):
    def __init__(self, num_features, embedding_dim=64):
        """
        num_features = dimension of each node's feature vector.
                       If each channel is a node, and you pass
                       [nChan, nTime], then num_features = nTime.
        """
        super(GraphAutoencoder, self).__init__()
        # Example with two ARMAConv layers
        self.encoder_gcn1 = pyg_nn.ARMAConv(num_features, 128, 
                                            num_stacks=2, 
                                            num_layers=3,
                                            shared_weights=True)
        self.encoder_gcn2 = pyg_nn.ARMAConv(128, embedding_dim, 
                                            num_stacks=2, 
                                            num_layers=3,
                                            shared_weights=True)
        self.decoder_fc1 = nn.Linear(embedding_dim, 128)
        self.decoder_fc2 = nn.Linear(128, num_features)

    def encode(self, x, edge_index):
        x = torch.relu(self.encoder_gcn1(x, edge_index))
        latent = torch.relu(self.encoder_gcn2(x, edge_index))
        return latent
    
    def decode(self, latent):
        x = torch.relu(self.decoder_fc1(latent))
        reconstructed = self.decoder_fc2(x)
        return reconstructed
    
    def forward(self, x, edge_index):
        latent = self.encode(x, edge_index)
        reconstructed = self.decode(latent)
        return reconstructed

# ===========================================
# 6) Training and evaluation
# ===========================================
def train_model_all(model, train_loader, device, num_epochs):
    model.to(device)
    optimizer = optim.Adam(model.parameters(), lr=1e-3)
    criterion = nn.SmoothL1Loss()

    for epoch in range(num_epochs):
        model.train()
        epoch_loss = 0.0
        for batch_data in train_loader:
            batch_data = batch_data.to(device)
            optimizer.zero_grad()

            # Forward
            reconstructed = model(batch_data.x, batch_data.edge_index)
            loss = criterion(reconstructed, batch_data.x)

            # Backprop
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()

        epoch_loss /= len(train_loader)
        print(f"Beta Model Training - Epoch [{epoch+1}/{num_epochs}], Loss: {epoch_loss:.4f}")
    return model

def evaluation_model(model, data_sample, edge_index, device):
    """
    data_sample: shape [nChan, nTime], or [nTime, nChan],
                 whichever you used in training.
    Returns: (original_np, reconstructed_np, error_np)
    """
    model.eval()
    with torch.no_grad():
        x = torch.tensor(data_sample, dtype=torch.float32, device=device)
        eidx = edge_index.to(device)
        reconstructed = model(x, eidx)
    orig = x.cpu().numpy()
    recon = reconstructed.cpu().numpy()
    error = orig - recon
    return orig, recon, error

# ===========================================
# 7) Main pipeline
#    1) For each beta, build dataset/loader.
#    2) Initialize & train a separate model.
#    3) Reconstruct each subject's data for that beta.
#    4) Store reconstruction in a final array.
# ===========================================

num_features = nTime   # if each channel is a node, and feature vector = [nTime]
num_epochs   = 200     # number of training - set as needed
batch_size   = 4      # dataset/batch to update hyperparameters
embedding_dim = 32  # latent dimension size

models = []
reconstructions_list = []  # will hold arrays of shape [nSubject, nChan, nTime] for each beta

if saveGAE == "xai":
    module_ratio = 0.3 # drop ratio used in module ablation
    drop_ratios = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] # drop ratios list used in random edge ablation
    ori_error_list = []; linPro_error_list = [] # linear probing
    imp_chan_list = []; imp_region_list = []; imp_time_list = [] # masking
    imp_latent_list = []; imp_module_list = []; imp_rand_edge_list = [] # ablation
    imp_node_edge_list = []; imp_chan_drop_list = [] #ablation
    am_list = []  # activation maximization
    

for iBeta in range(nBeta):
    print(f"\n=== Building dataset/model for Beta {iBeta} / {nBeta} ===")
    # Build dataset & loader
    full_dataset_i  = build_concat_dataset_for_beta(datain, iBeta)
    graph_dataset_i = EEGGraphDataset(full_dataset_i, edge_index)
    full_loader_i   = GeometricDataLoader(graph_dataset_i, batch_size=batch_size, shuffle=True)

    # Init model & train
    model_i = GraphAutoencoder(num_features=num_features, embedding_dim=embedding_dim)
    model_i = train_model_all(model_i, full_loader_i, device, num_epochs=num_epochs)
    models.append(model_i)
    
    if saveGAE == "xai":
        print(f"=== Explainable AI processes for Beta {iBeta} / {nBeta} ===")
        
        # Train linear probing layer
        linPro_model = xGAE.train_linear_probe(
            model_i, full_loader_i, device, embedding_dim=embedding_dim, num_features=num_features, num_epochs=num_epochs)
        
        # Apply Ablation 
        imp_latent = xGAE.apply_latent_ablation(model_i, full_loader_i, device, embedding_dim=embedding_dim) # latent ablation importance
        baseline_loss, imp_module = xGAE.apply_module_ablation(model_i, full_loader_i, device, xGAE.modules_names, ratio=module_ratio) # module ablation importance
        imp_rand_edge = xGAE.apply_random_edge_ablation(model_i, full_dataset_i, device, edge_index, drop_ratios, batch_size, baseline_loss) # random edge ablation importance
        imp_node_edge = xGAE.apply_node_edge_ablation(model_i, nChan, full_dataset_i, device, edge_index, batch_size, baseline_loss) # node-based edge ablation importance
        imp_latent_list.append(imp_latent); imp_module_list.append(imp_module)
        imp_rand_edge_list.append(imp_rand_edge); imp_node_edge_list.append(imp_node_edge)
        
        # Define helper lists
        ori_error_subjects = []; linPro_error_subjects = [] # original and linear probing errors
        imp_chan_subjects = []; imp_region_subjects = []; imp_time_subjects = [] # masking importance
        imp_chan_drop_subjects = [] # channel ablation importance
        am_list_latents = []  # activation maximization per latent dimension
        
        # Apply Activation Maximization for each latent dimension
        for latent_idx in range(embedding_dim):
            print(f"\n=== Activation Maximization for Beta {iBeta} / {nBeta}, Latent {latent_idx+1} ===")
            x_am = xGAE.apply_activation_maximization(
                model=model_i, edge_index=edge_index,
                nChan=nChan, nTime=nTime, device=device,
                latent_ind=latent_idx, steps=400,
                lr=0.01, l2=1e-3, tv=1e-3) # shape: [nChan, nTime]
            am_list_latents.append(np.array(x_am))
    
    # Reconstruct data for all subjects
    # each subject is an item in full_dataset_i
    # shape of full_dataset_i[s][0] = [nChan, nTime]
    recon_per_subject = [] 
    
    for s in range(len(full_dataset_i)):
        data_s = full_dataset_i[s][0].numpy()  # [nChan, nTime]
        _, recon_s, error_s = evaluation_model(model_i, data_s, edge_index, device)
        recon_per_subject.append(recon_s)  # shape [nChan, nTime]
        
        if saveGAE == "xai":
            # Original AE error
            ori_error_subjects.append(error_s)
            # Linear Probing error
            _, _, error_l = xGAE.evaluation_linear_probe(
                model_i, linPro_model, data_s, edge_index, device)
            linPro_error_subjects.append(error_l)
            
            # Masking
            imp_chan = xGAE.apply_channel_masking(model_i, data_s, edge_index, device, mode='mean') # mean method
            imp_region = xGAE.apply_region_masking(model_i, data_s, edge_index, device, mode='global_mean', k=1) # global mean method
            imp_time = xGAE.apply_time_frame_masking(model_i, data_s, edge_index, device, n_segments=100, mode='mean') # mean method, nSeg=100
            imp_chan_subjects.append(imp_chan); imp_region_subjects.append(imp_region); imp_time_subjects.append(imp_time)
            
            # Channel ablation
            imp_chan_drop = xGAE.apply_channel_dropping(model_i, data_s, edge_index, device)
            imp_chan_drop_subjects.append(imp_chan_drop)
    

    recon_per_subject = np.array(recon_per_subject)  # shape [nSubject, nChan, nTime]
    reconstructions_list.append(recon_per_subject)
    
    if saveGAE == "xai":
        ori_error_list.append(ori_error_subjects); linPro_error_list.append(linPro_error_subjects) # linear probing
        imp_chan_list.append(imp_chan_subjects); imp_region_list.append(imp_region_subjects); imp_time_list.append(imp_time_subjects) # masking
        imp_chan_drop_list.append(imp_chan_drop_subjects) # channel ablation
        am_list.append(np.array(am_list_latents))  # activation maximization
        print(f"=== Explainable AI processes are finished Beta {iBeta} / {nBeta} ===\n")

# Convert list -> array: shape [nBeta, nSubject, nChan, nTime]
reconstructed_all = np.array(reconstructions_list)

if saveGAE == "xai":
    ori_error_all = np.array(ori_error_list)           # shape [nBeta, nSubject, nChan, nTime]
    linPro_error_all = np.array(linPro_error_list)     # shape [nBeta, nSubject, nChan, nTime]
    imp_chan_all = np.array(imp_chan_list)             # shape [nBeta, nSubject, nChan]
    imp_region_all = np.array(imp_region_list)         # shape [nBeta, nSubject, nChan]
    imp_time_all = np.array(imp_time_list)             # shape [nBeta, nSubject, nSeg]
    imp_latent_all = np.array(imp_latent_list)         # shape [nBeta, embedding_dim]
    imp_rand_edge_all = np.array(imp_rand_edge_list)   # shape [nBeta, len(drop_ratios)]
    imp_node_edge_all = np.array(imp_node_edge_list)   # shape [nBeta, nChan]
    imp_chan_drop_all = np.array(imp_chan_drop_list)   # shape [nBeta, nSubject, nChan]
    am_all = np.array(am_list)                         # shape [nBeta, embedding_dim, nChan, nTime]


# ===========================================
# 8) Save XAI results if needed
# ===========================================

# Save reconstructed betas, original betas, adjacency matrix, model weights
if saveGAE == "yes" or saveGAE == "xai":
    savemat(os.path.join(newdir, "reconstructed_betas.mat"), {"reconstructed": reconstructed_all})
    savemat(os.path.join(newdir, "original_betas.mat"), {"data": datain})
    savemat(os.path.join(newdir, "adjacency_matrix.mat"), {"adjacency_matrix": binatry_matrix})
    for iBeta, model_i in enumerate(models):
        model_path = os.path.join(model_weights_dir, f"Model_Weights_Beta_{iBeta+1}.pt")
        torch.save(model_i.state_dict(), model_path) # Save trained model
        print(f"Model weights for Beta {iBeta+1} are saved.")
if saveGAE == "xai":
    savemat(os.path.join(LP_dir, "ori_error_all.mat"), {"ori_error_all": ori_error_all})
    savemat(os.path.join(LP_dir, "linPro_error_all.mat"), {"linPro_error_all": linPro_error_all})
    savemat(os.path.join(MK_dir, "imp_chan_all.mat"), {"imp_chan_all": imp_chan_all})
    savemat(os.path.join(MK_dir, "imp_region_all.mat"), {"imp_region_all": imp_region_all})
    savemat(os.path.join(MK_dir, "imp_time_all.mat"), {"imp_time_all": imp_time_all})
    savemat(os.path.join(AB_dir, "imp_latent_all.mat"), {"imp_latent_all": imp_latent_all})
    # use np.save for lists of dicts
    np.save(os.path.join(AB_dir, f'imp_module_list.npy'), np.array(imp_module_list))
    savemat(os.path.join(AB_dir, "imp_rand_edge_all.mat"), {"imp_rand_edge_all": imp_rand_edge_all})
    savemat(os.path.join(AB_dir, "imp_node_edge_all.mat"), {"imp_node_edge_all": imp_node_edge_all})
    savemat(os.path.join(AB_dir, "imp_chan_drop_all.mat"), {"imp_chan_drop_all": imp_chan_drop_all})
    savemat(os.path.join(AM_dir, "am_all.mat"), {"am_all": am_all})
    
# ===========================================
# 9) Viaualize XAI results if needed
# ===========================================
if saveGAE == "xai":
    xGAE.visualize_LP_result(ori_error_all, linPro_error_all, LP_dir) # linear probing results
    xGAE.visualize_CR_masking_importance(imp_chan_all, MK_dir, 'mean', target='channel') # channel masking results
    xGAE.visualize_CR_masking_importance(imp_region_all, MK_dir, 'global mean', target='region') # region masking results
    xGAE.visualize_time_importance(imp_time_all, MK_dir, 'mean', tmin_ms=-200, tmax_ms=600, tick_step_ms=50) # time frame masking results
    xGAE.visualize_module_ablation(imp_module_list, AB_dir, 0.3, xGAE.modules_names) # module ablation results
    xGAE.visualize_ablation_results(imp_latent_all, AB_dir, "Latent dimension", "Latent Ablation", xticklabels=None) # latent ablation results
    xGAE.visualize_ablation_results(imp_rand_edge_all, AB_dir, "Drop Ratio", "Random Edge Ablation", 
                                    xticklabels=[f"{int(ratio*100)}%" for ratio in drop_ratios]) # random edge ablation results
    xGAE.visualize_ablation_results(imp_node_edge_all, AB_dir, "Node (Channel)", "Node-based Edge Ablation", xticklabels=None) # node-based edge ablation results
    xGAE.visualize_ablation_results(imp_chan_drop_all, AB_dir, "Channel", "Channel Drop Ablation", xticklabels=None) # channel drop ablation results
    xGAE.plot_AM_beta_clusters_similarity(am_all, AM_dir, xGAE.beta_type) # activation maximization - beta PCA, cosine similarity and hierarchical dendrogram
    xGAE.plot_AM_channel_clusters(am_all, AM_dir) # activation maximization - channel clustering

# ===========================================
# 10) Final output to MATLAB
# ===========================================
learned_betas = {
    # shape [nBeta, nSubject, nChan, nTime]
    "reconstructed": reconstructed_all
}