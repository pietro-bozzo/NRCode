import torch.nn as nn


def loadModel(name,**kw_args):
    match name:
        case "CNN1d":
            model = CNN1D(**kw_args)
        case "CNN1D_noreg":
            model = CNN1D_noreg(**kw_args)
        case "CNN1d_dilation":
            model = CNN1D_dilation(**kw_args)
    return model


class CNN1D(nn.Module):
    def __init__(self,n_features):
        super().__init__()

        self.conv = nn.Sequential(
            nn.Conv1d(in_channels=n_features,out_channels=16,kernel_size=1), # PCA analogous
            nn.BatchNorm1d(16),
            nn.ReLU(),

            nn.Conv1d(in_channels=16,out_channels=32,kernel_size=31,padding=15),
            nn.BatchNorm1d(32),
            nn.ReLU(),

            nn.Conv1d(in_channels=32,out_channels=32,kernel_size=31,padding=15),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.MaxPool1d(4),

            nn.Conv1d(in_channels=32,out_channels=64,kernel_size=15,padding=7),
            nn.BatchNorm1d(64),
            nn.ReLU(),

            nn.Conv1d(in_channels=64,out_channels=64,kernel_size=15,padding=7),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.MaxPool1d(4),

            nn.AdaptiveAvgPool1d(128))

        self.fc = nn.Sequential(
            nn.Flatten(),
            nn.Linear(64*128,64),
            nn.ReLU(),
            nn.Linear(64,1),
            nn.Softplus()) # softplus for nonnegative target

    def forward(self,x):
        x = self.conv(x)
        x = x.squeeze(-1)
        return self.fc(x).squeeze(-1)


class CNN1D_noreg(nn.Module):
    def __init__(self,n_features):
        super().__init__()

        self.conv = nn.Sequential(
            nn.Conv1d(in_channels=n_features,out_channels=16,kernel_size=1), # PCA analogous
            nn.ReLU(),

            nn.Conv1d(in_channels=16,out_channels=32,kernel_size=31,padding=15),
            nn.ReLU(),

            nn.Conv1d(in_channels=32,out_channels=32,kernel_size=31,padding=15),
            nn.ReLU(),

            nn.Conv1d(in_channels=32,out_channels=64,kernel_size=15,padding=7),
            nn.ReLU(),

            nn.Conv1d(in_channels=64,out_channels=64,kernel_size=15,padding=7),
            nn.ReLU(),

            nn.AdaptiveAvgPool1d(16))

        self.fc = nn.Sequential(
            nn.Flatten(),
            nn.Linear(64*16,64),
            nn.ReLU(),
            nn.Linear(64,1),
            nn.Softplus()) # softplus for nonnegative target

    def forward(self,x):
        x = self.conv(x)
        x = x.squeeze(-1)
        return self.fc(x).squeeze(-1)


class CNN1D_dilation(nn.Module):
    def __init__(self,n_features):
        super().__init__()

        self.conv = nn.Sequential(
            nn.Conv1d(in_channels=n_features,out_channels=16,kernel_size=1), # PCA analogous
            nn.BatchNorm1d(16),
            nn.ReLU(),
            
            nn.Conv1d(in_channels=16,out_channels=32,kernel_size=15,padding=7,dilation=1),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            
            nn.Conv1d(in_channels=32,out_channels=32,kernel_size=15,padding=7,dilation=2),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.MaxPool1d(4),
            
            nn.Conv1d(in_channels=32,out_channels=64,kernel_size=15,padding=7,dilation=4),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            
            nn.Conv1d(in_channels=64,out_channels=64,kernel_size=15,padding=7,dilation=8),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.MaxPool1d(4),

            nn.AdaptiveAvgPool1d(16))

        self.fc = nn.Sequential(
            nn.Flatten(),
            nn.Linear(64*16,64),
            nn.ReLU(),
            nn.Linear(64,1),
            nn.Softplus()) # softplus for nonnegative target

    def forward(self,x):
        x = self.conv(x)
        x = x.squeeze(-1)
        return self.fc(x).squeeze(-1)