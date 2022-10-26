# HW2 问题记录

**c++ RBF神经网络**

1.  网络结构：
   1. 输入为一个x， 中间一个隐层激活函数为gauss函数，然后输出为一个y_pred
   2. loss为MSE
   3. optimal是沿梯度反向下降lr * grad
2. 参数
   1. W1：输入层到hidden layer的weight ,  shape = 1 * rbf_number
   2. B1: 输入层到hidden layer的bias, shape = 1 * rbf_number
   3. W2: hidden layer 到 输出y_pred的weight shape = 1 * rbf_number
   4. B2： 输出层的bias， shape = 1* 1
3.  forward:
   1. 公式：
      1. $Z = W^1 * x + B^1$
      2. $A = g_{0,1}(Z)$
      3. $y\_pred = \sum (A * W^2) + B^2$
   2. 其中向量使用ArrayXf进行**bitwise**加法和乘法
4.  backward：
   1.  公式：
      1.  $loss = (y - y\_pred)^2$
      2.  $DYpred = -(y - y\_pred)$
      3.  $DW^2 = DYpred * A$    //参数梯度
      4.  $DB = DYpred * 1$        //参数梯度
      5.  $DA = DYpred * W^2$
      6.  $DZ = DA * A * -Z$
      7.  $DW^1 = DZ * x$           //参数梯度
      8.  $DB^1 = DZ * 1$            //参数梯度
5.  问题描述：
   1.  使用imgui将采样的点作为输入进行训练，但是得到的参数进行绘制的图像是一条水平直线(输出是一个常数)