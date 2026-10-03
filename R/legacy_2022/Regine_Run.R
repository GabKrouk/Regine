###first make a Regine plot to see if your list of gene indeed share regulators more than expected by chance.
Regine.plot.gb(gene.list=i(), rep.per.range.value=100, kernel.range=2:9, name.your.list="Test.list.gab")


###check the pdf made in the folder and pick your kernel size accordingly

Find.me.direct.regulators.KERNEL.STAT.gb(gene.list=i(),kernel.size=6,rep=1000,  name.your.list="Test.regine.yoyoyo.dec.2022")

### 2 files containing your results will be made in the folder DAP.seq.live.report and DAP.seq.results.
