import matplotlib.pyplot as plt
import pandas as pd

waveform = pd.read_csv(r'./Inlet Flow Waveform/SV_70/Waveform_70_02.csv', names=['wave'])
print(waveform)
plt.plot(waveform['wave'])
plt.show()