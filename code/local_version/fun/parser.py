#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Project: Medulloblastoma Cavalli Kaggle Project
Description: Class to parse the data to local_version
Author: Thomas Neff
Date:  2024-11-01 11:57:33
"""

# import lib
#import pandas as pd
import modin.pandas as pd
import numpy as np
import os
from collections import OrderedDict

class Data():

    def __init__(self,exp_mat:pd.DataFrame=None,meta:pd.DataFrame=None,supp:OrderedDict=None):
        '''
        self.exp_mat: (pd.DataFrame) RNA expression matrix, genes x samples
        self.meta: (pd.DataFrame) metadata, samples x features
        self.supp: (OrderedDict)  supplementary data
        '''
        self.exp_mat=exp_mat or pd.DataFrame()
        self.meta=meta or pd.DataFrame()
        self.supp=supp or OrderedDict()

    def update_exp_mat(self,exp_mat:pd.DataFrame):
        '''
        Update the instance self.exp_mat

        exp_mat: (pd.DataFrame) RNA expression matrix, genes x samples
        '''
        self.exp_mat=exp_mat

    def update_meta(self,meta:pd.DataFrame):
        '''
        Update the instance self.meta

        meta: (pd.DataFrame) metadata, samples x features
        '''
        self.meta=meta

    def add_exp_mat(self, path_exp_mat:str, **kwargs):
        '''
        path_exp_mat: (str) RNA expression matrix pathway
        '''

        if not os.path.exists(path_exp_mat):
            raise ValueError(f"File {path_exp_mat} does not exist")

        # default values for kwargs
        kwargs.setdefault('sep','\t')
        kwargs.setdefault('on_bad_lines','warn')

        self.exp_mat=pd.read_csv(path_exp_mat,**kwargs)

    def add_meta(self, path_meta:str, **kwargs):
        '''
        Read the metadata and fill to self.meta instance

        path_meta: (str) metadata pathway
        '''

        if not os.path.exists(path_meta):
            raise ValueError(f"File {path_meta} does not exist")

        self.meta=pd.read_csv(path_meta,**kwargs)

    def add_genes_name(self,genes_name:list=[],exp_mat:pd.DataFrame=None,inplace:bool=True) -> pd.DataFrame:
        '''
        Add genes name annotation in RNA expression matrix (genes x samples)

        genes_name: (list) list of genes name annotation
        exp_mat: (pd.DataFrame) RNA expression matrix, genes x samples
        inplace: (bool) if inplace = True the instance self.exp_mat will be update otherwise the updated matrix will be returned, brings usability
        '''

        if self.exp_mat is None and exp_mat is None:
            raise ValueError(f"Provide RNA expression matrix (genes x samples)")

        if exp_mat is None:
            exp_mat = self.exp_mat

        if len(genes_name)!=exp_mat.shape[0] or len(genes_name)!=len(set(genes_name)):
            raise ValueError(f"The genes names list mismatching with the expression matrix rows")

        # add genes name list as rownames
        exp_mat.index=genes_name

        if inplace:
            self.update_exp_mat(exp_mat)
        else:
            return exp_mat

    def add_samples_name(self,samples_name:list=[],meta:pd.DataFrame=None,inplace:bool=True):
        '''
        Add samples name annotation in metadata (samples x fetaures)

        samples_name: (list) list of samples name annotation
        meta: (pd.DataFrame) metadata, samples x features
        inplace: (bool) if inplace = True the instance self.meta will be update otherwise the updated metadata will be returned
        '''

        if self.meta is None and meta is None:
            raise ValueError(f"Provide metadata (samples x features)")

        if meta is None:
            meta=self.meta

        if len(samples_name)!=meta.shape[0] or len(samples_name)!=len(set(samples_name)):
            raise ValueError(f"The samples names list mismatching with the  metadata rows")

        # add samples list as metadata rownames
        meta.index=samples_name

        if inplace:
            self.update_meta(meta)
        else:
            return meta


    @staticmethod
    def search_duplicates(values:list) -> set:
        '''
        Search for duplicate values in a list of values

        values (list) list of values
        '''

        return set([x for x in values if values.count(x) >1])

    def add_supp(self,key:str,value,inplace=True):
        '''
        Add supplementary data in OrderedDict

        key (str) key to OrderedDict
        value: (any) data/variables to store in object
        inplace (bool) replace the key args if the key is already in OrderedDict
        '''
        #  del set
        if key not in self.supp.keys() or inplace:
            self.supp[key]=value
        else:
            raise KeyError(f"The key {key} already exist in self.supp instance")

    def sel_samples_and_genes_name(self,samples_name:list=[],genes_name:list=[],exp_mat:pd.DataFrame=None,meta:pd.DataFrame=None,inplace:bool=True):
        '''
        Select the genes and also select the samples in metadata as well as expression matrix in same time

        samples_name: (list) list of samples name annotation
        genes_name: (list) list of genes name annotation
        exp_mat: (pd.DataFrame) RNA expression matrix, genes x samples
        meta: (pd.DataFrame) metadata, samples x features
        inplace: (bool) if inplace = True the instance self.meta and self.exp_mat will be update otherwise the updated metadata and expressiion matrix will be returned
        '''

        if self.exp_mat is None and exp_mat is None:
            raise ValueError(f"Provide RNA expression matrix (genes x samples)")

        if self.meta is None and meta is None:
            raise ValueError(f"Provide metadata (samples x features)")

        if meta is None:
            meta=self.meta

        if exp_mat is None:
            exp_mat=self.exp_mat

        # Select genes and samples name
        meta=meta.loc[samples_name,]

        if len(samples_name) > 0:
            exp_mat=exp_mat[samples_name]
        elif len(genes_name) > 0:
            exp_mat=exp_mat.loc[genes_name,]

        if inplace:
            self.meta=meta
            self.exp_mat=exp_mat
        else:
            return exp_at, meta
