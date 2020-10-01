
##### 剪切位点之后

start = 0
length = 50

func = get_region_icSHAPE_mean_relative_Splice
x, y, names = func(top80_isoform_list, start, length)
#j = sns.regplot(x, y)
j = sns.kdeplot(x, y, cmap="Greens", shade=True)
coor, p = scipy.stats.spearmanr(x, y)
Figures.annotate(plt.gca(), f"spearmanr={coor:.3} P={p:.3e}\nN={len(x)}",'top left')

x, y, names = func([d for d in canonical_isoform_list ], start, length) # if d.annot!='ORF7b'
for x_,y_,name in zip(x, y, names):
    plt.plot(x_, y_, '.', c=Colors.RGB['green'])
    plt.text(x_, y_, name)

plt.xlabel("icSHAPE reacvitity in 5’region of TRS-B")
plt.ylabel("log10(number of sgRNA)")
plt.title(f"start={start}, length={length}")
plt.show()


j = sns.regplot(x, y)
for x_,y_,name in zip(x, y, names):
    plt.text(x_, y_, name)

coor, p = scipy.stats.spearmanr(x, y)
j.annotate(s=f"spearmanr={coor:.3} P={p:.3e}\nN={len(x)}",xy=[ plt.xlim()[0]+0.05, plt.ylim()[1]-0.5 ])
plt.xlabel("icSHAPE reacvitity in 5’region of TRS-B")
plt.ylabel("log10(number of sgRNA)")
plt.title(f"start={start}, length={length}")
plt.savefig(join(HOME, "figs/1.pdf"))
plt.show()


##### 剪切位点之前

start = -40
length = 40

func = get_region_icSHAPE_mean_relative_Splice
x, y, names = func(top80_isoform_list, start, length)
#j = sns.regplot(x, y)
j = sns.kdeplot(x, y, cmap="Greens", shade=True)
coor, p = scipy.stats.spearmanr(x, y)
Figures.annotate(plt.gca(), f"spearmanr={coor:.3} P={p:.3e}\nN={len(x)}",'top left')

x, y, names = func([d for d in canonical_isoform_list if d.annot!='ORF7b'], start, length)
for x_,y_,name in zip(x, y, names):
    plt.plot(x_, y_, '.', c=Colors.RGB['green'])
    plt.text(x_, y_, name)

plt.xlabel("icSHAPE reacvitity in 5’region of TRS-L")
plt.ylabel("log10(number of sgRNA)")
plt.title(f"start={start}, length={length}")
plt.show()


j = sns.regplot(x, y)
for x_,y_,name in zip(x, y, names):
    plt.text(x_, y_, name)

coor, p = scipy.stats.spearmanr(x, y)
j.annotate(s=f"spearmanr={coor:.3} P={p:.3e}\nN={len(x)}",xy=[ plt.xlim()[0]+0.05, plt.ylim()[1]-0.5 ])
plt.xlabel("icSHAPE reacvitity in 5’region of TRS-L")
plt.ylabel("log10(number of sgRNA)")
plt.title(f"start={start}, length={length}")
plt.savefig(join(HOME, "figs/2.pdf"))
plt.show()

##### 与TE的关系

func = get_region_icSHAPE_mean_relative_Splice
x, y, names = func(top80_isoform_list, -50, 50, TE=TE_05h)
sns.regplot(x, y)
for x_,y_,name in zip(x, y, names):
    plt.text(x_, y_, name)

coor, p = scipy.stats.spearmanr(x, y)
Figures.annotate(plt.gca(), f"spearmanr={coor:.3} P={p:.3e}\nN={len(x)}",'top left')

plt.xlabel("icSHAPE reacvitity in flanking region of splice site")
plt.ylabel("TE (Finkel et al. RNA, 2020)")
plt.title("start=-50, length=50")
plt.show()





