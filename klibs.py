import pandas as pd, numpy as np
import pandas_profiling

print("!pip install dtale")
print("import pandas as pd, numpy as np")
def explore(df):
  '''
  explore(df) - uses pandas profiling to show lots of info about df
  '''
  display(pandas_profiling.ProfileReport(df, explorative=True))
