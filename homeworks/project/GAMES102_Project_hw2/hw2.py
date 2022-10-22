import torch
from torch import nn
from torch.utils.data import DataLoader, Dataset

import numpy as np
import matplotlib.pyplot as plt


lr = 0.001
epochs = 10
batch_size = 20


class RBFKernel(nn.Module):
    def __init__(self):
        super(RBFKernel, self).__init__()    
    def forward(self, x):
        return torch.exp(-0.5 * x * x) / (np.sqrt(2 * np.pi))
        
    
class RBFNetwork(nn.Module):
    def __init__(self, data_dim, layer_dim):
        super(RBFNetwork, self).__init__()
        self.data_dim = data_dim
        self.layer_dim = layer_dim
        
        #1. 先将x进行仿射变换(一层Linear)
        #2. 高斯激活函数， 构造高斯基函数空间
        #3. 使用基函数，线性组合成要训练的函数
        self.linear_gauss_stack = nn.Sequential(
            nn.Linear(self.data_dim, self.layer_dim),   
            RBFKernel(),
            nn.Linear(self.layer_dim, self.data_dim, bias=False)
        )
        self.linear = nn.Linear(data_dim, data_dim)
    def forward(self, x):
        res = self.linear_gauss_stack(x)
        #res = self.linear(x)
        return res

class RBF(nn.Module):
    def __init__(self, rbf_number):
        super(RBF, self).__init__()
        #自定义神经网络参数a, b (要训练的两组参数)
        self.rbf_number = rbf_number
        self.activate_func = RBFKernel()
        
        self.linear = nn.Linear(rbf_number, 1, bias = True)
        
        
        self.linear_weight = nn.Parameter(torch.ones(rbf_number), requires_grad=True)
        self.linear_bias = nn.Parameter(torch.ones(1), requires_grad=True)
        
        
        #参数矩阵转置，到row sapce进行计算
        self.a = nn.Parameter(torch.ones(rbf_number),requires_grad=True)
        self.b = nn.Parameter(torch.ones(rbf_number), requires_grad=True)
        self.init()
        
    def init(self): 
        self.a.data.normal_(0, 0.2)
        self.b.data.normal_(0, 0.2)
        self.linear_weight.data.normal_(0, 0.2)
        
        #self.linear.weight.data.normal_(0, 0.2)        
        
    def forward(self, x):
        p = self.activate_func(self.a * x + self.b)
        
        #方法1 bit-wise 乘法，加和等价与方法2中的matmult,这两种方法就是参数矩阵的维度以及size设置不同
        #self.linear_weight = nn.Parameter(torch.ones(rbf_number), requires_grad=True)
        #([p1, p2, p3, p4] * [w1, w2, w3,w4]).sum() + b === 向量相乘
        y = (self.linear_weight * p).sum() + self.linear_bias
        
        #方法2
        #self.linear_weight = nn.Parameter(torch.ones(1, rbf_number), requires_grad=True)
        #y = torch.matmul(self.linear_weight, p) + self.linear_bias
        
        #方法3 library中的linear
        #y = self.linear(p)
        #print(y)
        return y
        
def training(data, model, loss_fn, optimizer):
    
    for epoch in range(epochs):
        print(f"epoch: {epoch}")
        #training
        for x, y in data.items():
            xt = torch.Tensor([x])
            yt = torch.Tensor([y])
            
            #预测，计算loss
            y_predict = model(xt)
            
            loss = loss_fn(y_predict, yt)
            
            #优化
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            # print(f'a.grad = {model.a.grad}\n b.grad={model.b.grad}')
            # for name, parameters in model.named_parameters():
            #     print(f"{name}, {parameters}")
            # print("============================")

def testing(model, test_data):
    #testing
    with torch.no_grad():
        for x, y in test_data.items():
            xt = torch.Tensor([x])
            yt = torch.Tensor([y])
            y_pred = model(xt)
            print(f'{xt}->{y_pred}->{yt}')
    

class PointDataSet(Dataset):
    def __init__(self, sample_data, label, batch_size):
        self.sample_data = sample_data
        self.label = label
        self.batch_size = batch_size
    def __len__(self):
        return len(self.sample_data)
    def __getitem__(self, idx):
        return self.sample_data[idx], self.label[idx]


if __name__ == "__main__":
    
    #model = RBFNetwork(1, 1000)
    model = RBF(10) 
    loss_fn = nn.MSELoss().to('cuda')
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    
    
    for name, parameters in model.named_parameters():
        print(f"{name}, {parameters}")
    
    data = {x:y for x, y in zip(range(100), range(100))}
    
    training(data, model, loss_fn, optimizer)
    
    #用model，绘图
    data_num= 100
    xs = np.linspace(0, 10, data_num)
    ys = []
    with torch.no_grad():
        for x in xs:
        
            xt = torch.Tensor([x])
            y = model(xt).sum().item()
            ys.append(y)
    #plot
    fig, ax = plt.subplots()  

    ax.plot(xs, ys, linewidth=1.0)

    plt.show()      
    
    # x = torch.tensor([1, 2, 3])
    # y = torch.tensor([1,2,3])   

    # data_set = PointDataSet(x, y, 1)
    
    # data_loader = DataLoader(data_set)
    # model = RBFNetwork(1, 3)
    # loss_fn = nn.MSELoss()
    # optimizer = torch.optim.SGD(model.parameters(), lr=lr)
    
    # for x, y in enumerate(data_loader):
    #     #预测,计算loss
    #     y_predict = model(x)    
    #     loss = loss_fn(y_predict, y)
        
    #     #求梯度
    #     optimizer.zero_grad()
    #     loss.backward()
        
    #     #梯度下降
    #     optimizer.setp()
        
    
    # training(data_loader, model, loss_fn, optimizer)
 

    
                
        
        