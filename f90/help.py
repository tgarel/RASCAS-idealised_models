from scipy.io import FortranFile as ff

f=ff('help','r')
flux = f.read_reals()
f.close()


print(flux.size)
print(flux[flux==0].size)
