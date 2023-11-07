# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 01:45:00 2023

@author: Tingpeng Yang
"""

import torch
import torch.nn as nn
from .transform import EmbeddingTransform

class LogisticActivation(nn.Module):
    """
    Implementation of Generalized Sigmoid
    Applies the element-wise function:

    :math:`\\sigma(x) = \\frac{1}{1 + \\exp(-k(x-x_0))}`

    :param x0: The value of the sigmoid midpoint
    :type x0: float
    :param k: The slope of the sigmoid - trainable -  :math:`k \\geq 0`
    :type k: float
    :param train: Whether :math:`k` is a trainable parameter
    :type train: bool
    """

    def __init__(self, x0=0, k=1, train=False):
        super(LogisticActivation, self).__init__()
        self.x0 = x0
        if train:
            self.k = nn.Parameter(torch.FloatTensor([float(k)]))
        else:
            self.k = k

    def forward(self, x):
        x = torch.clamp(1 / (1 + torch.exp(-self.k * (x - self.x0))), min=0, max=1).squeeze()
        return x
        
class ModelInteraction(nn.Module):
    def __init__(
        self,
        args,
        contact,
        do_sigmoid=True,
        gamma_init=0,
        p0=0.5
    ):
        super(ModelInteraction, self).__init__()
        self.one_word = args.one_word
        self.device=args.device
        self.do_sigmoid = do_sigmoid
        if do_sigmoid:
            self.activation = LogisticActivation(x0=p0, k=20, train=True)
        self.embedding_transform_s = EmbeddingTransform(
            args.input_dim, args.projection_dim, words=5,dropout=args.dropout_p,
            nhead=args.nhead,num_layers=args.num_layers,
        )
        self.embedding_transform_l = EmbeddingTransform(
            args.input_dim, args.projection_dim, words=4**args.one_word + 1,dropout=args.dropout_p,
            nhead=args.nhead,num_layers=args.num_layers,
        )
        self.contact = contact
        self.gamma = nn.Parameter(torch.FloatTensor([gamma_init]))

    def map_predict(self, z0, z1):
        z1=self.embedding_transform_s(z1)
        z0=self.embedding_transform_l(z0)
        C = self.contact.forward(z0, z1)
        yhat = C
        # Mean of contact predictions where p_ij > mu + gamma*sigma
        mu = torch.mean(yhat,dim=(1,2,3)).repeat(yhat.shape[2]*yhat.shape[3],1).T.reshape(yhat.shape[0],yhat.shape[1],yhat.shape[2],yhat.shape[3])
        sigma = torch.var(yhat,dim=(1,2,3)).repeat(yhat.shape[2]*yhat.shape[3],1).T.reshape(yhat.shape[0],yhat.shape[1],yhat.shape[2],yhat.shape[3])
        Q = torch.relu(yhat - mu - (self.gamma * sigma))
        phat = torch.sum(Q,dim=(1,2,3)) / (torch.sum(torch.sign(Q),dim=(1,2,3)) + 1)
        if self.do_sigmoid:
            phat = self.activation(phat).reshape(phat.shape[0])
        return C, phat