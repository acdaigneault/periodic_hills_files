import os


ext = '.csv'

for i in range(1, 10):
    os.rename(r'C:\Users\Acdai\OneDrive - polymtl.ca\Polytechnique\Session A2020\Données\Reynolds5600\Rapp2009\Rapp2009_' + str(i) + ext,
              r'C:\Users\Acdai\OneDrive - polymtl.ca\Polytechnique\Session A2020\Données\Reynolds5600\Rapp2009\Rapp2009_0' + str(i) + ext)


