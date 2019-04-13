import pandas as pd
import Smith_Waterman

DNA=pd.read_csv("DNA.csv")
print(DNA)

Smith_Waterman.test()
