#!/usr/bin/env python3
# NiPyAEoutliers.py

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, ConcatDataset
from torch_geometric.loader import DataLoader as GeometricDataLoader
from torch_geometric.data import Data
from torch_geometric.utils import dense_to_sparse
import torch_geometric.nn as pyg_nn

# ===========================================
# 1) Retrieve inputs from MATLAB
# ===========================================
#  MATLAB calls:
#    learned_betas = pyrunfile("NiPyAEoutliers.py", "learned_betas",
#                              datain=beta_values, binatry_matrix=neighbourgh_matrix);
#

# ===========================================
datain = np.array(datain)          # shape: [nBeta, nTime, nChan, nSubject] (adjust as needed)
binatry_matrix = np.array(binatry_matrix)

# Parse shapes
nBeta     = datain.shape[0]
nTime     = datain.shape[1]
nChan     = datain.shape[2]
nSubject  = datain.shape[3]



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

    def forward(self, x, edge_index):
        x = torch.relu(self.encoder_gcn1(x, edge_index))
        latent = torch.relu(self.encoder_gcn2(x, edge_index))
        x = torch.relu(self.decoder_fc1(latent))
        reconstructed = self.decoder_fc2(x)
        return reconstructed

# ===========================================
# 6) Training and evaluation
# ===========================================
def train_model_all(model, train_loader, device, num_epochs=300):
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

num_features = nTime  # if each channel is a node, and feature vector = [nTime]
num_epochs   = 300     # set as needed

models = []
reconstructions_list = []  # will hold arrays of shape [nSubject, nChan, nTime] for each beta

for iBeta in range(nBeta):
    print(f"\n=== Building dataset/model for Beta {iBeta} ===")
    # Build dataset & loader
    full_dataset_i = build_concat_dataset_for_beta(datain, iBeta)
    graph_dataset_i = EEGGraphDataset(full_dataset_i, edge_index)
    full_loader_i = GeometricDataLoader(graph_dataset_i, batch_size=4, shuffle=True)

    # Init model & train
    model_i = GraphAutoencoder(num_features=num_features, embedding_dim=32)
    model_i = train_model_all(model_i, full_loader_i, device, num_epochs=num_epochs)
    models.append(model_i)

    # Reconstruct data for all subjects
    # each subject is an item in full_dataset_i
    # shape of full_dataset_i[s][0] = [nChan, nTime]
    recon_per_subject = []
    for s in range(len(full_dataset_i)):
        data_s = full_dataset_i[s][0].numpy()  # [nChan, nTime]
        _, recon_s, _ = evaluation_model(model_i, data_s, edge_index, device)
        recon_per_subject.append(recon_s)  # shape [nChan, nTime]

    recon_per_subject = np.array(recon_per_subject)  # shape [nSubject, nChan, nTime]
    reconstructions_list.append(recon_per_subject)

# Convert list -> array: shape [nBeta, nSubject, nChan, nTime]
reconstructed_all = np.array(reconstructions_list)

# ===========================================
# 8) Final output to MATLAB
# ===========================================
learned_betas = {
    # shape [nBeta, nSubject, nChan, nTime]
    "reconstructed": reconstructed_all
}