import pandas as pd, numpy as np
import dtale
print("!pip install dtale")
print("import pandas as pd, numpy as np")
def explore(df):
  '''
  explore(df) - uses dtale to show lots of info about df
  '''
  dtale.show(df)
