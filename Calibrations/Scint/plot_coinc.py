#!/usr/bin/python

#######################################################
# 5/1/20
# John Williamson
#
# Plots coincidence peak widths from different calibrations
#
#######################################################


import argparse
#from DB_reader import DB_reader
import pandas
import matplotlib.pyplot as plt

from matplotlib.ticker import ScalarFormatter


# define arguments for script

parser = argparse.ArgumentParser(description='plot_coinc.py: Plots coindcidence peak widths from different calibrations (from chosen csv file).')
parser.add_argument("--csv_name", default="Defaulto", help="Name of csv file to be read in by script")


args = parser.parse_args()


csv_name = args.csv_name


coinc_df = pandas.read_csv(csv_name)

print(coinc_df)

print("\n")




print(f"{coinc_df['Description'].tolist()}")

print("\n")
print(f"{type(coinc_df['Description'].tolist())}")



fig, axs = plt.subplots()

plt.figure(dpi=1200)

Desc_list = coinc_df['Description'].tolist()

width_list = coinc_df['Width'].tolist()

# y_formatter = ScalarFormatter(useOffset=False)
# axs.yaxis.set_major_formatter(y_formatter)


width_list, Desc_list = (list(t) for t in zip(*sorted(zip(width_list, Desc_list))))


print(f"width_list type = {type(width_list)}")
print(f"width_list[0] = {type(width_list[0])}")

width_list = [t*(1e9) for t in width_list]
#width_list = width_list*(1e9)


axs.set_ylabel(r'Coincidence Peak, $\sigma$ [ns]')
axs.set_xlabel(r'Calibration method')


axs.grid(True)
axs.set_axisbelow(True)

#axs.plot(Desc_list,width_list,'.', color = ['Black','Red','Blue'])
axs.scatter(Desc_list,width_list,marker='.', s= 200, color=['Red','Blue','Green','Orange'])





#plt.show()



#fig.savefig(f"ME_plots/{plot_order}_order/{elem['Element']}_{fname_ext}.png", bbox_inches='tight', pad_inches=0)

fig.savefig(f"coinc_csv/plots/coinc_timing.pdf")
