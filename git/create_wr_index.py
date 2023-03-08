import pandas as pd

#Converting weather regimes from numbers to names and vice versa
def era_indexer(era_type):
  main_dir = '../../'
  if era_type=='era5':
    input = "../../data/bundle/regimes/Z0500_N161_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local.txt" #snakemake.input.wr
    converter = {0:"No", 1:"AT", 2:"AR", 3:"GL", 4:"EuBL", 5:"ScBL", 6:"ZO", 7:"ScTr"}
  #elif era_type=='erai':
  #  input = '%s/data_ERAI/bundle/regimes/Z0500_N81_Atl_EU2_year_6h_7_10_7_ncl_all_LCO_local.txt'%(main_dir)
  #  converter = {0:"No", 1:"AT", 2:"GL", 3:"AR", 4:"ScTr", 5:"EuBL", 6:"ScBL", 7:"ZO"}
  input_file = pd.read_csv(input, skiprows=7, header=None, delim_whitespace=True, index_col=[1])[4].squeeze()
  input_file.index = pd.to_datetime(input_file.index, format='%Y%m%d_%H')
  output_file = pd.Series([converter[x] for x in input_file], index=input_file.index)
  return(output_file)

def era_reindexer(input):
  era5_reconvert = {"No":0, "AT":1, "AR":2, "GL":3, "EuBL":4, "ScBL":5, "ZO":6, "ScTr":7}
  output_file = pd.Series([era5_reconvert[x] for x in input], index=input.index)
  return(output_file)

if __name__ == "__main__":
    era5 = era_indexer('era5')
    era5.to_csv("../../data/bundle/wr5_index.txt") #snakemake.output.wr_index)