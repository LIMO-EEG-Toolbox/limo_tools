#!/usr/bin/env python3
# Explainable_GAE.py
# Yanhaoyang Li

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import os
import matplotlib.pyplot as plt
from torch_geometric.loader import DataLoader as GeometricDataLoader
from torch_geometric.utils import k_hop_subgraph, subgraph
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import cosine_similarity
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib.colors import BoundaryNorm
from collections import defaultdict
from torch_geometric.data import Data

# ===========================================
# 1) Helping Functions and Variables
# ===========================================

# The map from channel index to channel name
channel_name = {
    0:  "FP1", 1:  "F3", 2:  "F7", 3:  "FC3", 4:  "C3", 5:  "C5",
    6:  "P3", 7:  "P7", 8:  "P9", 9:  "PO7", 10: "PO3", 11: "O1",
    12: "Oz", 13: "Pz", 14: "CPz", 15: "FP2", 16: "Fz", 17: "F4",
    18: "F8", 19: "FC4", 20: "FCz", 21: "Cz", 22: "C4", 23: "C6", 
    24: "P4", 25: "P8", 26: "P10", 27: "PO8", 28: "PO4", 29: "O2"
}

# List of module names for module ablation
modules_names = [
        "encoder_gcn1",
        "encoder_gcn2",
        "decoder_fc1",
        "decoder_fc2",
    ]

# The map from beta index to experimental setting
# The number is the nBeta in each paradigm
beta_type = {2: {0 : "Standard", 1: "Deviant"}, # MMN
             4: {0: "Cars", 1:"Faces" , 2 : "Cars-scrambled", 3: "Faces-scrambled"}, # N170
             8: {0: "prime,related,list1", 1: "prime,related,list2", 2: "prime,unrelated,list1",
                      3: "prime,unrelated,list2", 4: "target,related,list1", 5: "target,related,list2", 
                      6: "target,unrelated,list1", 7: "target,unrelated,list2"}, # N400
             25: {0: "AA target", 1: "AB non-target", 2: "AC non-target", 3: "AD non-target", 4: "AE non-target", 
                    5: "BA non-target", 6: "BB target", 7: "BC non-target", 8: "BD non-target", 9: "BE non-target",
                    10: "CA non-target", 11: "CB non-target", 12: "CC target", 13: "CD non-target", 14: "CE non-target", 15: "DA non-target", 
                    16: "DB non-target", 17: "DC non-target", 18: "DD target", 19: "DE non-target", 
                    20: "EA non-target", 21: "EB non-target", 22: "EC non-target", 23: "ED non-target", 24: "EE target"} # P3
             }

# Compute k-hop neighborhood channels for each channel
def compute_neighborhood_channels(edge_index, k=1):
    neighborhood_channels = []
    for c in range(edge_index.max().item() + 1):
        # Get k-hop subgraph for channel c
        channels, _, _, _ = k_hop_subgraph(c, k, edge_index)
        neighborhood_channels.append(channels.cpu().numpy())
    return neighborhood_channels

# Compute baseline loss for ablation studies
def compute_loss(model, data_loader, device):
    model.eval()
    total_loss, n = 0, 0
    criterion = nn.SmoothL1Loss()
    with torch.no_grad():
        # Iterate over data loader
        for data in data_loader:
            data = data.to(device)
            recon = model(data.x, data.edge_index)
            loss = criterion(recon, data.x)
            total_loss += loss.item()
            n += 1
    return total_loss / max(n, 1)

# Define EEG Graph Dataset for PyG
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
# 2) Linear Probing
# ===========================================

# Define a simple linear layer probe without activation function as decoder
class LinearProbe(nn.Module):
    def __init__(self, embedding_dim, num_features):
        super().__init__()
        # Linear layer to map from embedding_dim to num_features
        self.probe = nn.Linear(embedding_dim, num_features, bias=True)

    def forward(self, x):
        return self.probe(x)


# Train linear probe using frozen GAE latent representations
def train_linear_probe(model, train_loader, device, embedding_dim, num_features, num_epochs):
    # Freeze GAE model parameters
    model.eval()
    for param in model.parameters():
        param.requires_grad = False
    
    # Initialize linear probe
    linear_probe = LinearProbe(embedding_dim, num_features).to(device)
    optimizer = optim.Adam(linear_probe.parameters(), lr=1e-3)
    criterion = nn.SmoothL1Loss()
    
    for epoch in range(num_epochs):
        linear_probe.train()
        epoch_loss = 0.0
        for batch_data in train_loader:
            batch_data = batch_data.to(device)
            optimizer.zero_grad()
            
            # Get latent representations from frozen GAE
            with torch.no_grad():
                latent = model.encode(batch_data.x, batch_data.edge_index)
            
            lin_reconstructed = linear_probe(latent)
            loss = criterion(lin_reconstructed, batch_data.x)
             
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()

        epoch_loss /= len(train_loader)
        print(f"Linear Probe Training - Epoch [{epoch+1}/{num_epochs}], Loss: {epoch_loss:.4f}")

    return linear_probe

def evaluation_linear_probe(model, linear_probe, data_sample, edge_index, device):
    model.eval()
    linear_probe.eval()
    with torch.no_grad():
        x = torch.tensor(data_sample, dtype=torch.float32, device=device)
        eidx = edge_index.to(device)
        # Get latent representations
        latent_repr = model.encode(x, eidx)
        # Reconstruct using linear probe
        lin_reconstructed = linear_probe(latent_repr)
    orig = x.cpu().numpy()
    recon = lin_reconstructed.cpu().numpy()
    error = orig - recon
    return orig, recon, error

# ===========================================
# 3) Masking
# ===========================================
#   A. Channel Masking
#   B. Channel-based Region Masking
#   C. Time Frame Masking
# ===========================================

###### Channel Masking ######

# Define channel masking function with 'zero' and 'mean' methods
def channel_masking(x_n, ch_idx, mode='mean'):
    x_masked = x_n.copy()
    if mode == 'zero':
        x_masked[ch_idx, :] = 0.0
    # Mean of all other channels
    elif mode == 'mean':
        nChan = x_n.shape[0]
        mask = np.ones(nChan, dtype=bool)
        mask[ch_idx] = False
        channel_mean = x_n[mask, :].mean(axis=0)
        x_masked[ch_idx, :] = channel_mean
    return x_masked

# Apply channel masking and compute importance
def apply_channel_masking(model, sample_n, edge_index, device, mode='mean'):
    model.eval()
    with torch.no_grad():
        x = torch.tensor(sample_n, dtype=torch.float32, device=device)
        recon = model(x, edge_index.to(device))
        # Compute baseline loss
        base_loss = nn.SmoothL1Loss()(recon, x).item()
        
        nChan = sample_n.shape[0]
        imp = np.zeros(nChan)
        
        # Iterate over channels to compute importance
        # The higher the loss increase, the more important the channel
        for c in range(nChan):
            x_masked = channel_masking(sample_n, c, mode=mode)
            x_masked = torch.tensor(x_masked, dtype=torch.float32, device=device)
            recon_masked = model(x_masked, edge_index.to(device))
            # Masked loss is computed between recon_masked and original x
            masked_loss = nn.SmoothL1Loss()(recon_masked, x).item()
            imp[c] = masked_loss - base_loss
            
    return imp

###### Channel-based region masking ######

# Define region masking function with 'zero', 'global_mean', and 'local_mean' methods
def region_masking(x_n, neig_list, mode='local_mean'):
    x_masked = x_n.copy()

    # Get the rest of the channels except the neighborhood channels and itself
    nChan = x_n.shape[0]
    mask_neig = np.zeros(nChan, dtype=bool)
    mask_neig[neig_list] = True
    mask_rest = ~mask_neig
    
    if mode == 'zero':
        x_masked[mask_neig, :] = 0.0
    # The mean of the rest of the channels
    elif mode == 'global_mean':
        if mask_rest.any():
            rest_mean = x_n[mask_rest, :].mean(axis=0)
        # If no rest channels, use overall mean
        else:
            rest_mean = x_n.mean(axis=0)
        x_masked[mask_neig, :] = rest_mean
    # The mean of the neighborhood channels and itself
    elif mode == 'local_mean':
        if mask_rest.any():
            local_mean = x_n[mask_neig, :].mean(axis=0)
        else:
            local_mean = x_n.mean(axis=0)
        x_masked[mask_neig, :] = local_mean
    return x_masked

# Apply region masking and compute importance
def apply_region_masking(model, sample_n, edge_index, device, mode='mean', k=1):
    model.eval()
    with torch.no_grad():
        # Precompute neighborhood channels for all channels
        neighborhood_channels = compute_neighborhood_channels(edge_index, k)
        x = torch.tensor(sample_n, dtype=torch.float32, device=device)
        recon = model(x, edge_index.to(device))
        base_loss = nn.SmoothL1Loss()(recon, x).item()

        nChan = sample_n.shape[0]
        imp = np.zeros(nChan)
        
        for c in range(nChan):
            # Get neighborhood channels list for channel c
            neig_list = neighborhood_channels[c]
            x_masked = region_masking(sample_n, neig_list, mode=mode)
            x_masked = torch.tensor(x_masked, dtype=torch.float32, device=device)
            recon_masked = model(x_masked, edge_index.to(device))
            masked_loss = nn.SmoothL1Loss()(recon_masked, x).item()
            imp[c] = masked_loss - base_loss
    return imp

###### Time frame masking ######

# Define time frame masking function with 'zero' and 'mean' methods
def time_frame_masking(x_n, seg_idx, boundaries, mode='mean'):
    x_masked = x_n.copy()
    nTime = x_n.shape[1]
    
    # Get start and end indices for the segment
    start = boundaries[seg_idx]
    end = boundaries[seg_idx + 1]
    # Set the values in the time frame to zero or mean
    if mode == 'zero':
        x_masked[:, start:end] = 0.0
    # The mean of all other time points
    elif mode == 'mean':
        mask = np.ones(nTime, dtype=bool)
        mask[start:end] = False
        time_mean = x_n[:, mask].mean(axis=1, keepdims=True)
        x_masked[:, start:end] = time_mean
    return x_masked

# Apply time frame masking and compute importance
def apply_time_frame_masking(model, sample_n, edge_index, device, n_segments=100, mode='mean'):
    model.eval()
    with torch.no_grad():
        x = torch.tensor(sample_n, dtype=torch.float32, device=device)
        recon = model(x, edge_index.to(device))
        base_loss = nn.SmoothL1Loss()(recon, x).item()
        # Get the boundaries for time segments
        nTime = sample_n.shape[1]
        imp = np.zeros(n_segments)
        boundaries = np.linspace(0, nTime, n_segments + 1, dtype=int)
        
        for seg_idx in range(n_segments):
            x_masked = time_frame_masking(sample_n, seg_idx, boundaries, mode=mode)
            x_masked = torch.tensor(x_masked, dtype=torch.float32, device=device)
            recon_masked = model(x_masked, edge_index.to(device))
            masked_loss = nn.SmoothL1Loss()(recon_masked, x).item()
            imp[seg_idx] = masked_loss - base_loss
    return imp
        
            
# ===========================================
# 4) Ablation
# ===========================================
#  A. Latent Ablation
#  B. Module Ablation
#  C. Edge Ablation
#  D. Channel Ablation
# ===========================================

###### Latent ablation Part ######

def apply_latent_ablation(model, data_loader, device, embedding_dim):
    model.eval()
    with torch.no_grad():
        # Compute baseline loss
        baseline_loss = compute_loss(model, data_loader, device)
        imp = np.zeros(embedding_dim)
        criterion = nn.SmoothL1Loss()
        
        for dim in range(embedding_dim):
            print(f"Ablating latent dimension {dim+1}/{embedding_dim}")
            total_loss, n = 0, 0
            # Iterate over data loader
            for data in data_loader:
                data = data.to(device)
                latent = model.encode(data.x, data.edge_index)
                latent_ablated = latent.clone()
                # Ablate one dimension using zeroing method
                latent_ablated[:, dim] = 0  
                recon = model.decode(latent_ablated)
                loss = criterion(recon, data.x)
                total_loss += loss.item()
                n += 1
            # Compute average loss over all batches
            avg_loss = total_loss / max(n, 1)
            imp[dim] = avg_loss - baseline_loss
    return imp

###### Module ablation ######

def compute_loss_module_ablation(model, data_loader, device, module_name, ratio=0.3):
    # Get the target module from the model
    modules = dict(model.named_modules())
    if module_name not in modules:
        raise ValueError(f"Module {module_name} is not found in model.")
    target_module = modules[module_name]
    mask = None
    
    # Define forward hook to ablate the output of the target module
    def hook(module, input, output):
        # Make mask only once
        nonlocal mask
        
        if mask is None:
            # Get number of units in the output
            num = output.size(-1)
            # Determine number of units to ablate
            k = int(num * ratio)
            # Randomly select k units to ablate
            ind = torch.randperm(num, device=output.device)[:k]
            mask = torch.ones(num, device=output.device)
            mask[ind] = 0.0
        return output * mask
    
    # Register the forward hook
    handle = target_module.register_forward_hook(hook)
    ablated_loss = compute_loss(model, data_loader, device)
    handle.remove()
    return ablated_loss

# Apply module ablation for a list of modules and compute importance
def apply_module_ablation(model, full_loader_i, device, modules_name, ratio=0.3):
    baseline_loss = compute_loss(model, full_loader_i, device)
    module_ablation_res = {}
    # Iterate over modules to compute ablation importance
    for mod_name in modules_name:
        ablated_loss = compute_loss_module_ablation(model, full_loader_i, device, mod_name, ratio=ratio)
        module_ablation_res[mod_name] = ablated_loss - baseline_loss
    return baseline_loss, module_ablation_res

###### Edge ablation ######

# Randomly remove a ratio of edges from the graph and compute loss increase
def make_random_edge_ablation(edge_index, drop_ratio):
    edge_index = edge_index.cpu()
    # Get number of edges
    edge_num = edge_index.size(1)
    # Determine number of edges to drop
    num_drop = int(edge_num * drop_ratio)
    # Randomly select edges to drop
    ind = torch.randperm(edge_num)
    mask = torch.ones(edge_num, dtype=torch.bool)
    mask[ind[:num_drop]] = False
    edge_index_dropped = edge_index[:, mask]
    return edge_index_dropped

# Remove edges conneted to a selected node
def make_node_edge_ablation(edge_index, node_idx):
    edge_index = edge_index.cpu()
    # Create mask to exclude edges connected to the node
    mask = (edge_index[0] != node_idx) & (edge_index[1] != node_idx)
    edge_index_dropped = edge_index[:, mask]
    return edge_index_dropped

# Apply random edge ablation and compute importance
def apply_random_edge_ablation(model, full_dataset_i, device, edge_index, drop_ratios, batch_size, baseline_loss):
    imp_rand_edge = []
    for dr in drop_ratios:
        # Generate randomly ablated edge index
        edge_index_rand = make_random_edge_ablation(edge_index, dr)
        # Create new dataset and data loader with ablated edges index
        rand_graph_dataset_i = EEGGraphDataset(full_dataset_i, edge_index_rand)
        rand_full_loader_i = GeometricDataLoader(rand_graph_dataset_i, batch_size=batch_size, shuffle=True)
        edge_rand_loss = compute_loss(model, rand_full_loader_i, device)
        # Relative importance
        imp_rand_edge.append((edge_rand_loss - baseline_loss) / (baseline_loss+1e-12))
    return np.array(imp_rand_edge)

# Apply node-based edge ablation and compute importance
def apply_node_edge_ablation(model, nChan, full_dataset_i, device, edge_index, batch_size, baseline_loss):
    imp_node_edge = []
    for c in range(nChan):
        # Generate node-ablated edge index
        edge_index_node_ablated = make_node_edge_ablation(edge_index, c)
        # Create new dataset and data loader with ablated edges index
        node_graph_dataset_i = EEGGraphDataset(full_dataset_i, edge_index_node_ablated)
        node_full_loader_i = GeometricDataLoader(node_graph_dataset_i, batch_size=batch_size, shuffle=True)
        edge_node_loss = compute_loss(model, node_full_loader_i, device)
        imp_node_edge.append(edge_node_loss - baseline_loss)
    return np.array(imp_node_edge)


###### Channel Ablation ######

# Define channel dropping method
def apply_channel_dropping(model, sample_n, edge_index, device):
    model.eval()
    with torch.no_grad():
        x = torch.tensor(sample_n, dtype=torch.float32, device=device)
        e_idx = edge_index.to(device)
        recon = model(x, e_idx)

        nChan = sample_n.shape[0]
        imp = np.zeros(nChan)
        all_channels = torch.arange(nChan, device=device)
        
        for c in range(nChan):
            mask = torch.ones(nChan, dtype=torch.bool, device=device)
            mask[c] = False
            masked_channels = all_channels[mask]
            # Create subgraph with masked channels, dropping the channel c and its edges from the graph
            base_loss_kept = nn.SmoothL1Loss()(recon[masked_channels, :], x[masked_channels, :]).item()
            masked_edge_index = subgraph(masked_channels, e_idx, relabel_nodes=True)[0]
            # Drop the channel c from the data x
            masked_x = x[masked_channels, :]
            recon_masked = model(masked_x, masked_edge_index)
            masked_loss = nn.SmoothL1Loss()(recon_masked, masked_x).item()
            imp[c] = masked_loss - base_loss_kept
    return imp

# ===========================================
# 5) Activation Maximization
# ===========================================

# Define activation maximization function
def activation_maximization(
    activation_func, input_shape, device,
    steps=400, lr=0.01, l2=1e-3, tv=1e-3):
    
    # Initialize input with random noise
    x = torch.randn(input_shape, device=device, requires_grad=True)
    optimizer = torch.optim.Adam([x], lr=lr)
    
    # Optimization loop
    for step in range(steps):
        optimizer.zero_grad()
        # Compute activation
        activation = activation_func(x)
        
        # Compute loss = negative activation + regularization
        loss = -activation
        loss += l2 * torch.mean(x**2)
        
        # Total Variation (TV) regularization
        dx = x[1:, :] - x[:-1, :]
        dy = x[:, 1:] - x[:, :-1]
        tv_loss = dx.abs().mean() + dy.abs().mean()
        loss += tv * tv_loss

        loss.backward()
        optimizer.step()
    
        if (step+1) % 50 == 0:
            print(f"Step {step+1} | loss {loss.item():.4f} | activation {activation.item():.4f}")
    
    return x.detach().cpu()

# Apply activation maximization to a specific latent dimension or overall latent space
def apply_activation_maximization(
    model, edge_index, nChan, nTime,
    device, latent_ind, steps=400, lr=0.01,
    l2=1e-3, tv=1e-3):
    
    model = model.to(device)
    model.eval()
    eidx = edge_index.to(device)
    
    # Define activation function
    def activation_func(x):
        latent = model.encode(x, eidx)
        # Activation is mean of specific latent dimension or overall latent space
        if latent_ind == None:
            act = latent.mean()
        else:
            act = latent[:, latent_ind].mean()
        return act
    
    # Perform activation maximization
    x = activation_maximization(
        activation_func=activation_func,
        input_shape=(nChan, nTime),
        device=device, steps=steps,
        lr=lr, l2=l2, tv=tv)
    
    return x

# ===========================================
# 6) Visualization Functions
# ===========================================
#  A. Visualize Linear Probing Results
#  B. Visualize channel and region masking importance
#  C. Visualize time frame masking importance
#  D. Visualize module ablation importance
#  E. Visualize latent dimension clustering and similarity based on AM results
#  F. Visualize channel clustering based on AM results
# ===========================================

# Visualize Linear Probing Results comparison with original GAE across time
def visualize_LP_result(ori_error_all, linPro_error_all, save_dir):
    ori_error_all = ori_error_all[:-1]
    linPro_error_all = linPro_error_all[:-1]
    # Determine global max and min for consistent y-axis limits
    global_max = np.percentile(np.abs(np.concatenate((ori_error_all, linPro_error_all))), 95)
    global_min = 0
    nBeta = ori_error_all.shape[0]

    # Print MAE of original GAE and Linear Probing and plot reconstruction 
    # result comparison of origainal GAE and Linear Probing
    for iBeta in range(nBeta):
        err_o = ori_error_all[iBeta]  # [nSubject, nChan, nTime]
        err_l = linPro_error_all[iBeta]
        # mean_err_o = np.mean(np.abs(err_o))
        # mean_err_l = np.mean(np.abs(err_l))
        # print(f"Beta {iBeta+1} Mean Absolute Error:")
        # print(f"Original GAE: {mean_err_o:.4f}")
        # print(f"Linear Probing: {mean_err_l:.4f}")
        
        time_mae_o = np.abs(err_o).mean(axis=(0, 1))
        time_mae_l = np.abs(err_l).mean(axis=(0, 1))
        
        plt.figure(figsize=(16,13))
        plt.plot(time_mae_o, label='Original GAE', linewidth=2)
        plt.plot(time_mae_l, label='Linear Probing', linewidth=2)
        plt.xlabel('Time')
        plt.ylabel('Error (MAE)')
        plt.yticks(np.arange(global_min, global_max, step=0.01))
        plt.ylim(global_min, global_max)
        plt.title(f'Beta {iBeta+1} error (MAE) over time')
        plt.legend(fontsize=40, loc='upper left')
        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, f'Beta_{iBeta+1}_LP.png'))
        plt.close()

# Visualize aggregated importance across Betas for masking
def visualize_masking_importance_across_betas(
    data,
    savepath,
    target,
    lablename='Channel Index',
    x_positions=None,      # bar positions
    x_ticks=None,          # tick positions
    x_ticklabels=None     # tick labels
):
    nBeta, nUnits = data.shape
    norm_all_imp = np.zeros_like(data, dtype=float)

    # Normalize per Beta
    for iBeta in range(nBeta):
        value = data[iBeta]
        max_abs = np.max(np.abs(value))
        norm_all_imp[iBeta] = value / max_abs if max_abs > 0 else value

    # Cross Beta Average
    agg_imp = norm_all_imp.mean(axis=0)

    # Bar positions
    if x_positions is None:
        x_positions = np.arange(nUnits)

    fig, ax = plt.subplots(figsize=(10, 6))

    # Bar width for time-like continuous x
    if lablename =='Time (ms)' and len(x_positions) > 1:
        width = (x_positions[1] - x_positions[0]) * 0.8
        ax.bar(x_positions, agg_imp, width=width, align='center')
    else:
        ax.bar(x_positions, agg_imp)

    if x_ticks is None:
        # default ticks match positions
        x_ticks = x_positions
        
    if x_ticklabels is None:
        # default labels match ticks
        x_ticklabels = [str(int(x)+1) for x in x_ticks]
        
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_ticklabels, fontsize=12)
    ax.tick_params(axis='x', labelsize=12, rotation=0)
    ax.tick_params(axis='y', labelsize=12)
    ax.set_xlabel(lablename, fontsize=16)
    ax.set_ylabel('Aggregated Importance (normalized)', fontsize=16)
    ax.set_title(f'{target.capitalize()} Masking Aggregated Importance across Betas',
                 fontsize=18, pad=10)

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', linestyle='--', alpha=0.4)

    fig.tight_layout()
    fig.savefig(savepath, dpi=300)
    plt.close(fig)


# Visualize importance of channel and region masking
def visualize_CR_masking_importance(data, save_path, method, target='channel'):
    """
    data: [nBeta, nSubject, nChan]
    target: 'channel' or 'region'
    """
    data = data[:-1] # exclude the last Beta
    nBeta, _, nChan = data.shape
    global_max = np.percentile(data, 97)
    global_min = np.min(data)
    # [nBeta, nChan], mean across subjects
    values = data.mean(axis=1) 
    
    visualize_masking_importance_across_betas(
        values,
        os.path.join(save_path, f'Aggregated_Importance_{target.capitalize()}_Masking_{method.capitalize()}.png'),
        target,
        lablename='Channel Index')

    # Heatmap across Betas and Channels
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(
        values,
        aspect='auto',
        origin='lower',
        cmap='Reds',
        vmin=global_min,
        vmax=global_max
    )

    # Set axis labels, font sizes, and title
    ax.set_xlabel(target.capitalize(), fontsize=16)
    ax.set_ylabel('Beta Index', fontsize=16)
    ax.set_title(f'{target.capitalize()} Masking - {method.capitalize()}',
                fontsize=18, pad=10)

    # Customize tick labels
    ax.set_xticks(np.arange(nChan))
    ax.set_xticklabels(np.arange(1, nChan + 1), fontsize=10)
    ax.set_yticks(np.arange(nBeta))
    ax.set_yticklabels(np.arange(1, nBeta + 1), fontsize=12)

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label('Importance', fontsize=14)
    cbar.ax.tick_params(labelsize=12)

    fig.tight_layout()
    fig.savefig(
        os.path.join(save_path, f'Heatmap_{target.capitalize()}_Masking_Importance_{method.capitalize()}.png'),
        dpi=300
    )
    plt.close(fig)


    for iBeta in range(nBeta):
        # Mean across subjects
        values = data[iBeta].mean(axis=0)

        fig, ax = plt.subplots(figsize=(10, 6))

        x = np.arange(nChan)
        ax.bar(x, values)

        # Set axis labels, font sizes, and title
        ax.set_xticks(x)
        ax.set_xticklabels(np.arange(1, nChan + 1))
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)

        ax.set_xlabel(f'Channel index', fontsize=18)
        ax.set_ylabel('Importance', fontsize=18)
        ax.set_ylim(global_min, global_max)
        ax.set_title(
            f'Beta {iBeta+1} - {target.capitalize()} Masking - {method.capitalize()}',
            fontsize=14,
            pad=12
        )

        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # Add grid lines for y-axis
        ax.grid(axis='y', linestyle='--', alpha=0.4)

        # Save figure
        fig.tight_layout()
        fig.savefig(os.path.join(save_path, f'Beta_{iBeta+1}_{target.capitalize()}_Masking_Importance_{method.capitalize()}.png'), dpi=300)
        plt.close(fig)
        
# Visualize time frame masking importance
def visualize_time_importance(
    data, 
    save_path, 
    method, 
    tmin_ms=-200, 
    tmax_ms=600, 
    tick_step_ms=50
    ):
    """
    data: [nBeta, nSubject, nSeg]
    x-axis shown in ms from tmin_ms to tmax_ms with ticks every tick_step_ms
    """
    data = data[:-1]
    nBeta, _, nSeg = data.shape
    global_max = np.percentile(data, 97)
    global_min = np.min(data)

    # [nBeta, nSeg] mean across subjects
    values = data.mean(axis=1)

    # build time segments in ms
    boundaries_ms = np.linspace(tmin_ms, tmax_ms, nSeg + 1)          # [nSeg+1]
    centers_ms = (boundaries_ms[:-1] + boundaries_ms[1:]) / 2.0      # [nSeg]
    seg_width_ms = (tmax_ms - tmin_ms) / nSeg
    xticks_ms = np.arange(tmin_ms, tmax_ms + 1e-9, tick_step_ms)

    # Heatmap across Betas and Time Segments
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(
        values,
        aspect='auto',
        origin='lower',
        cmap='Reds',
        vmin=global_min,
        vmax=global_max,
        extent=[tmin_ms, tmax_ms, 0, nBeta]   # make x-axis in ms
    )

    ax.set_xlabel('Time (ms)', fontsize=16)
    ax.set_ylabel('Beta Index', fontsize=16)
    ax.set_title(f'Time frame masking - {method.capitalize()}',
                 fontsize=18, pad=10)
    
    # Customize tick labels
    ax.set_xticks(xticks_ms)
    ax.tick_params(axis='x', labelsize=10)
    ax.set_yticks(np.arange(nBeta) + 0.5)
    ax.set_yticklabels(np.arange(1, nBeta + 1), fontsize=12)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label('Importance', fontsize=14)
    cbar.ax.tick_params(labelsize=12)

    fig.tight_layout()
    fig.savefig(
        os.path.join(save_path, f'Time_Frame_Masking_Importance_{method.capitalize()}.png'),
        dpi=300
    )
    plt.close(fig)

    xticks_ms = np.arange(-200, 601, 50)

    visualize_masking_importance_across_betas(
        values,
        os.path.join(save_path, f'Aggregated_Importance_Time_Frame_{method.capitalize()}.png'),
        lablename='Time (ms)',
        target='Time Frame',
        x_positions=centers_ms,
        x_ticks=xticks_ms,
        x_ticklabels=[str(int(x)) for x in xticks_ms]
    )

    # Plot per Beta importance
    for iBeta in range(nBeta):
        values_beta = data[iBeta].mean(axis=0)  # [nSeg]

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.bar(centers_ms, values_beta, width=seg_width_ms * 0.8)

        ax.set_xlim(tmin_ms, tmax_ms)
        ax.set_xticks(xticks_ms)
        ax.set_xticklabels([str(int(x)) for x in xticks_ms])

        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)

        ax.set_xlabel('Time (ms)', fontsize=18)
        ax.set_ylabel('Importance', fontsize=18)
        ax.set_ylim(global_min, global_max)
        ax.set_title(
            f'Beta {iBeta+1} - Time Frame Masking - {method.capitalize()}',
            fontsize=14,
            pad=12
        )

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', linestyle='--', alpha=0.4)

        fig.tight_layout()
        fig.savefig(
            os.path.join(save_path, f'Beta_{iBeta+1}_Time_Frame_Masking_Importance_{method.capitalize()}.png'),
            dpi=300
        )
        plt.close(fig)
        
# Visualize aggregated importance across Betas for module ablation
def visualize_ablation_importance_across_betas(
    data, 
    save_dir, 
    modules_names, 
    labelname,
    title,
    x_ticklabels=None):
    
    if labelname == 'Module Name':
        nBeta = len(data)
        nMod = len(modules_names)

        # Build matrix [nBeta, nMod]
        M = np.zeros((nBeta, nMod), dtype=float)
        for iBeta, d in enumerate(data):
            row = np.array([float(d.get(m, 0.0)) for m in modules_names], dtype=float)

            max_abs = np.max(np.abs(row))
            if max_abs > 0:
                row = row / max_abs
            M[iBeta, :] = row
    else:
        if data.ndim == 3:
            data = data.mean(axis=1)  # mean across subjects
        nBeta, nUnits = data.shape
        M = data.copy()
        # Normalize per Beta
        for iBeta in range(nBeta):
            row = M[iBeta]
            max_abs = np.max(np.abs(row))
            if max_abs > 0:
                M[iBeta] = row / max_abs
    
    if x_ticklabels is None:
        x_ticklabels = list(range(1, nUnits + 1))
        
        
    agg = M.mean(axis=0)
    x = np.arange(len(agg))
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(x, agg)
    
    ax.set_xticks(x)
    ax.set_xticklabels(x_ticklabels, fontsize=12)
    ax.tick_params(axis='y', labelsize=12)
    ax.set_xlabel(labelname, fontsize=16)
    ax.set_ylabel('Aggregated Importance (normalized)', fontsize=16)
    ax.set_title(f'{title} Aggregated Importance', fontsize=18, pad=10)
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.axhline(0, linewidth=1)

    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, f'Aggregated_Importance_{title}.png'), dpi=300)
    plt.close(fig)
        

# Visualize module ablation importance
def visualize_module_ablation(data, save_dir, module_ratio, modules_names):
    imp_module_list = data[:-1]
    nBeta = len(imp_module_list)

    global_max = np.max([max(imp.values()) for imp in imp_module_list])
    global_min = np.min([min(imp.values()) for imp in imp_module_list])

    visualize_ablation_importance_across_betas(
        imp_module_list, 
        save_dir, 
        modules_names, 
        labelname='Module Name',
        title='Module Ablation',
        x_ticklabels = modules_names
    )
    
    for iBeta in range(nBeta):
        module_imp_dict = imp_module_list[iBeta]

        modules = list(module_imp_dict.keys())
        values = np.array(list(module_imp_dict.values()))
        x = np.arange(len(values))
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.bar(x, values)

        # Set axis labels & ticks
        ax.set_xticks(x)
        ax.set_xticklabels(modules, fontsize=16)
        ax.tick_params(axis='y', labelsize=12)
        ax.set_xlabel('Module Name', fontsize=18)
        ax.set_ylabel('Importance', fontsize=16)
        ax.set_ylim(global_min, global_max)
        ax.set_title(f'Beta {iBeta+1}-Module importance (ratio = {module_ratio})', fontsize=18, pad=10)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', linestyle='--', alpha=0.4)
        ax.axhline(0, linewidth=1)

        fig.tight_layout()
        fig.savefig(
            os.path.join(
                save_dir,
                f'Beta_{iBeta+1}_Module_Ablation_Ratio_{module_ratio}.png'
            ),
            dpi=300
        )
        plt.close(fig)
    
# Visualize ablation results for latent, edge, and channel ablation        
def visualize_ablation_results(data, save_dir, xlabel, title, xticklabels=None):
    '''
    Shape of Latent Ablation Data:[nBeta, embedding_dim]
    Shape of node-based Edge Ablation Data: [nBeta, nChan]
    Shape of random Edge Ablation Data: [nBeta, len(drop_ratios)]
    Shape of Channel Ablation Data: [nBeta, nSubject, nChan]
    '''
    data = data[:-1]  # exclude the last Beta
    nBeta = data.shape[0]
    global_max = np.percentile(data, 97)
    global_min = np.percentile(data, 3)
    
    visualize_ablation_importance_across_betas(
        data, 
        save_dir, 
        modules_names=None, 
        labelname=xlabel,
        title=title,
        x_ticklabels=xticklabels
    )
    
    for iBeta in range(nBeta):
        if data.ndim == 3:
            value = data[iBeta].mean(axis=0)  # mean across subjects
        else:
            value = data[iBeta]
            
        n_units = len(value)
        x = np.arange(n_units)
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.bar(x, value)
        
        ax.set_xticks(np.arange(n_units))
        if xticklabels is None:
            ax.set_xticklabels(np.arange(1, n_units + 1), fontsize=12)
        else:
            ax.set_xticklabels(xticklabels, fontsize=12)
            
        ax.tick_params(axis='y', labelsize=12)
        ax.set_xlabel(xlabel, fontsize=16)
        ax.set_ylabel("Importance", fontsize=16)
        ax.set_ylim(global_min, global_max)
        ax.set_title(f'Beta {iBeta+1} - {title}', fontsize=18, pad=10)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', linestyle='--', alpha=0.4)
        ax.axhline(0, linewidth=1)
        
        fig.tight_layout()
        fig.savefig(os.path.join(save_dir, f'Beta_{iBeta+1}_{title}.png'), dpi=300)
        plt.close(fig)
        

# Visualize clustering of latent dimensions across betas based on AM results
def plot_AM_beta_clusters_similarity(data, save_dir, beta_type):
    # Prepare data for clustering
    all_x_am_use = data[:-1]
    nBeta, embedding_dim, nChan, nTime = all_x_am_use.shape
    units = all_x_am_use.reshape(nBeta*embedding_dim, nChan*nTime)  # shape: [nBeta*embedding_dim, nChan*nTime]

    # label ids for each unit
    beta_ids = np.repeat(np.arange(nBeta), embedding_dim)  # shape: [nBeta*embedding_dim]

    # Standardize data
    scaler = StandardScaler()
    units_std = scaler.fit_transform(units)

    # perform PCA to reduce dimensionality
    pca = PCA(n_components=10)
    units_pca = pca.fit_transform(units_std)

    # PCA for visualization
    pca2 = PCA(n_components=2)
    units_2d = pca2.fit_transform(units_pca)

    # Plot clustering results
    plt.figure(figsize=(7, 6))
    cmap = plt.get_cmap('tab10', nBeta)
    boundaries = np.arange(0.5, nBeta + 1.5, 1)
    norm = BoundaryNorm(boundaries, cmap.N)

    scatter = plt.scatter(units_2d[:, 0], units_2d[:, 1], c=beta_ids+1, cmap=cmap, norm=norm, alpha=0.6, s=25)

    beta_map = beta_type.get(nBeta, {})

    tick_values = np.arange(1, nBeta + 1)
    beta_indices = tick_values - 1
    tick_labels = [beta_map.get(int(b), str(int(b)+1)) for b in beta_indices]

    cbar = plt.colorbar(scatter)
    cbar.set_label("Beta")
    cbar.set_ticks(tick_values)
    cbar.set_ticklabels(tick_labels)

    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("PCA visualization of latent dimensions colored by beta")
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f'PCA_visualization_latent_dimensions.png'))
    plt.close() 
    
    beta_mean = all_x_am_use.mean(axis=1)
    x_beta = beta_mean.reshape(nBeta, nChan*nTime)  # shape: [nBeta, nChan*nTime]
    x_beta_std = scaler.fit_transform(x_beta)
    # Compute cosine similarity
    cos_sim = cosine_similarity(x_beta_std)

    beta_map = beta_type.get(nBeta, {})
    beta_labels = [beta_map.get(i, str(i+1)) for i in range(nBeta)]

    # Plot heatmap of cosine similarity
    plt.figure(figsize=(7, 6))
    x_sim = plt.imshow(cos_sim, vmin=-1, vmax=1)
    plt.colorbar(x_sim, label='Cosine similarity')

    plt.xticks(np.arange(nBeta), beta_labels, rotation=45, ha='right')
    plt.yticks(np.arange(nBeta), beta_labels)

    plt.xlabel('Beta')
    plt.ylabel('Beta')
    plt.title('Similarity between betas (mean AM)')
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f'Beta_similarity_heatmap.png'))


    # plot hierarchical clustering dendrogram of betas
    link = linkage(x_beta_std, method='ward')
    plt.figure(figsize=(7, 6))
    dendrogram(link, labels=beta_labels)
    plt.ylabel('Distance')
    plt.title('Hierarchical Dendrogram of Betas')
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f'Beta_hierarchical_dendrogram.png'))
    plt.close()

# Visualize clustering of channels based on AM results 
def plot_AM_channel_clusters(data, save_dir):
    data = data[:-1]
    # reshpe for channel clustering
    nChan = data.shape[2]
    x_chan = data.transpose(2, 0, 1, 3).reshape(nChan, -1)

    # Standardize data
    scaler = StandardScaler()
    x_chan_std = scaler.fit_transform(x_chan)

    # KMeans clustering
    k = 5  # according to the brain regions
    km = KMeans(n_clusters=k, random_state=0, n_init='auto')
    chan_labels = km.fit_predict(x_chan_std)

    pca = PCA(n_components=2)
    x_chan_2d = pca.fit_transform(x_chan_std)

    cmap = plt.get_cmap('tab10', k)
    boundaries = np.arange(0.5, k + 1.5, 1)
    norm = BoundaryNorm(boundaries, cmap.N)

    plt.figure(figsize=(7, 6))
    scatter = plt.scatter(x_chan_2d[:, 0], x_chan_2d[:, 1], c=chan_labels, cmap=cmap, norm=norm, alpha=0.6, s=25)

    cbar = plt.colorbar(scatter)
    cbar.set_label("Cluster")
    cbar.set_ticks(np.arange(1, k+1, 1))
    cbar.set_ticklabels(np.arange(1, k+1, 1))

    for ch in range(nChan):
        name = channel_name[ch]
        plt.text(x_chan_2d[ch, 0], x_chan_2d[ch, 1], name, fontsize=9)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('KMeans Clustering of Channels Based on AM')
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f'KMeans_clustering_channels.png'))
    plt.close()


    # group channel indices by cluster
    clusters = defaultdict(list)
    for ch_idx, cl in enumerate(chan_labels):
        clusters[cl].append(ch_idx)

    # print channels in each cluster
    for cl in range(k):
        print(f"\n=== Cluster {cl} ===")
        ch_list = clusters[cl]
        ch_list = sorted(ch_list)
        for ch_idx in ch_list:
            name = channel_name[ch_idx]
            print(f"  Channel {ch_idx+1:2d}  ->  {name}")

    

        