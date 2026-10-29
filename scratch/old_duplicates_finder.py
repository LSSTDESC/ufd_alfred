# old duplicate dropping code (it adds some things 2x)

moresig_Peaks = []
for i in range(len(Peaks)):
    overlaps = []
    for j in range(len(Peaks)):
        if i==j:
            continue
        else:
            angsep = projector.angsep(Peaks[i].ra, Peaks[i].dec, Peaks[j].ra, Peaks[j].dec)
            if angsep<Peaks[i].r:
                overlaps.append(Peaks[j])
    try:
        mostsig_Peak_i = max(overlaps, key=lambda x: x.sig)
        moresig_Peaks.append(mostsig_Peak_i)
    except:
        moresig_Peaks.append(Peaks[i])
#in case of duplicates (probably could make this logic better)
Peaks = list(set(moresig_Peaks))
print('Number of peaks after: ', len(Peaks))      
