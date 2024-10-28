#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Project: Medulloblastoma Cavalli Kaggle Project
Description: Class RLEplot
Author: Thomas Neff
Date: 2024-10-28 16:36:11
"""

# ADD COMMENTS

class RLEplot():

    def __init__(self,rle_mat:pd.DataFrame=None):

        self.rle_mat=rle_mat
        self.sample_stats=OrderedDict()

    @classmethod
    def add_rle_mat(cls,exp_mat:pd.DataFrame) -> "RLEplot":

        # calculate the median for each gene accross samples and calculate RLE values
        gene_medians=exp_mat.median(axis=1)
        rle_mat=exp_mat.sub(gene_medians,axis=0)

        return cls(rle_mat=rle_mat)

    def comp_summary_stats(self,rle_mat:pd.DataFrame=None,outliers:bool=True):

        if self.rle_mat is None and rle_mat is None:
            raise ValueError(f"Provide RLE matrix")

        if rle_mat is None:
            rle_mat=self.rle_mat

        sample_medians=[] # medians for each samples
        sample_q1=[] # first quartile for each samples
        sample_q3=[] # third quartile for each samples
        outlier_points={} # outliers for each samples
        samples=rle_mat.columns

        for sample in samples:
            rle_sample=rle_mat[sample]
            median=rle_sample.median()
            q1=rle_sample.quantile(0.25)
            q3=rle_sample.quantile(0.75)

            # append sample statistics
            sample_medians.append(median)
            sample_q1.append(q1)
            sample_q3.append(q3)

            # print(q1)
            # print(q3)

            # identify outliers based on IQR criterion
            if outliers:

                iqr=q3-q1
                outlier_values = rle_sample[(rle_sample < q1 - 1.5 * iqr) | (rle_sample > q3 + 1.5 * iqr)]
                outlier_points[sample]=outlier_values.to_list()

        self.sample_stats=self._fill_sample_stats(median=sample_medians,q1=sample_q1,q3=sample_q3,outliers=outlier_points,samples=samples)

    @staticmethod
    def _fill_sample_stats(median:list,q1:list,q3:list,outliers:list,samples:list) -> OrderedDict:
        return OrderedDict([('median',median),('q1',q1),('q3',q3),('outliers',outliers),('samples',samples)])

    def get_outliers_stats(self,rle_mat:pd.DataFrame=None) -> pd.DataFrame:

        if self.rle_mat is None and rle_mat is None:
            raise ValueError(f"Provide RLE matrix")

        if rle_mat is None:
            rle_mat=self.rle_mat

        if 'outliers' not in set(self.sample_stats.keys()):
            self.comp_summary_stats(rle_mat=rle_mat,outliers=True)

        samples=self.sample_stats['samples']
        n_out=[len(self.sample_stats['outliers'][sample])for sample in samples]
        n_genes=rle_mat.shape[0]
        p_out=[n/n_genes for n in n_out]

        return pd.DataFrame(data=np.array([n_out,p_out]),index=['# out','% out'],columns=samples)

    def get_summary_stats(self,rle_mat:pd.DataFrame=None,outliers:bool=True) -> pd.DataFrame:

        if self.rle_mat is None and rle_mat is None:
            raise ValueError(f"Provide RLE matrix")

        if rle_mat is None:
            rle_mat=self.rle_mat

        self.comp_summary_stats(rle_mat,outliers)

        index=['q1','median','q3']
        summary_stats=pd.DataFrame(data=np.array([self.sample_stats[i] for i in index]),index=index,columns=self.sample_stats['samples'])

        if outliers:
            return pd.concat([summary_stats,self.get_outliers_stats(rle_mat)])
        else:
            return summary_stats

    def density(self,rle_mat:pd.DataFrame=None,ref:bool=True,outliers:bool=False,rout:float=0.1,figsize=(14, 7),save:str=""):

        if self.rle_mat is None and rle_mat is None:
            raise ValueError(f"Provide RLE matrix")

        if rle_mat is None:
            rle_mat=self.rle_mat

        if rout <= 0 or rout > 1:
            raise ValueError(f'rout = {rout}, rout is a ratio and must be between 0 and 1')

        self.comp_summary_stats(rle_mat,outliers)

        # create plot
        plt.figure(figsize=figsize)

        # convert samples name to integer position
        sample_pos=range(len(self.sample_stats['samples']))

        # median
        plt.plot(sample_pos,self.sample_stats['median'],color="dimgrey",label="median")

        # add line for quartiles Q1 and Q3
        plt.plot(sample_pos, self.sample_stats['q1'], color="dodgerblue", label="Q1")
        plt.plot(sample_pos, self.sample_stats['q3'], color="crimson", label="Q3")

        # fill the area between Q1 and Q3
        plt.fill_between(sample_pos, self.sample_stats['q1'], self.sample_stats['q3'], color="lightgrey", alpha=0.5, label="IQR") # add median

        # add reference line at zero
        if ref:
            plt.plot(sample_pos, [0 for _ in range(len(sample_pos))], color="snow", linestyle="--", linewidth=1)

        if outliers and len(self.sample_stats['outliers'].keys())>0:

            for pos,sample in zip(sample_pos,self.sample_stats['samples']):
                outlier_values=self.sample_stats['outliers'][sample]
                # density is too high for an understanding representation
                n_outliers=len(outlier_values)
                for i in np.linspace(0,n_outliers-1,math.floor(n_outliers*rout)):
                    plt.plot(pos,outlier_values[math.floor(i)],'.', markersize=5, alpha=0.3, color='black')

        plt.xlabel("Samples")
        plt.ylabel("Relative Log Expression")
        plt.legend(loc="upper right")
        plt.show()

        if len(save) != 0:
            plt.savefig(save)
