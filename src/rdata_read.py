import pyreadr

result = pyreadr.read_r('/Users/danakeydar/IDC/Methylation-master/ch3/train_data/BRCA.RData') # also works for Rds

# done! let's see what we got
# result is a dictionary where keys are the name of objects and the values python
# objects
print(result.keys()) # let's check what objects we got
# df1 = result["df1"] # extract the pandas data frame for object df1

methyl = result['methyl']
expressi = result['expressi']

print(f" methyl shape: {methyl.shape}")
print(f" expressi shape: {expressi.shape}")

print(methyl.head(20))
print(expressi.head(20))