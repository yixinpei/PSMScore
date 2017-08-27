def getMassMap(massFileName):
	massFile = open(massFileName)
	massMap = {}
	for line in massFile:
		if not line=='\n':
			elems = line.split('=')
			aa = elems[0].strip()
			m = float(elems[1].strip())
			massMap[aa] = m
	massFile.close()
	return massMap

def getModifyMap(modifyFileName):
	modifyFile = open(modifyFileName)
	modifyMap = {}
	i = 0
	skip = 2
	for line in modifyFile:
		i = i+1
		if i>skip and not line=='\n' and line.find('->')==-1 and line.find('name')==-1:
			elems = line.split('=')
			name = elems[0]
			delMass = float(elems[1].split(' ')[2])
			modifyMap[name] = delMass
	modifyFile.close()
	return modifyMap

def getPeptideMQ(peptide, charge, modifyIndexes, modifyNames, massMap, modifyMap):
	masses = [massMap[p] for p in peptide]
	delmasses = [modifyMap[name] for name in modifyNames]
	for index, delmass in zip(modifyIndexes, delmasses):
		masses[index-1] += delmass
	curMasses=[masses[0]]
	for mass in masses[1:]:
		curMasses.append(curMasses[-1] + mass)
	b = [[mass/(q+1)+massMap['H'] for mass in curMasses] for q in xrange(charge)]
	reCurMasses = [masses[-1]]
	for mass in reversed(masses[:-1]):
        	reCurMasses.append(reCurMasses[-1] + mass)
	y = [[(mass+massMap['H2O'])/(q+1)+massMap['H'] for mass in reversed(reCurMasses)] for q in xrange(charge)]
	return b, y

def getSpectraData(spectraFileName):
	spectraFile = open(spectraFileName)
	mqs = []
	intensitys = []
	for line in spectraFile:
		if not line == '\n':
			elems = line.split(' ')
			mqs.append(float(elems[0]))
			intensitys.append(float(elems[1]))
	maxIntensity = max(intensitys)
	intensitys = [intensity/maxIntensity for intensity in intensitys]
	return mqs, intensitys
	spectraFile.close()
	return mqs, intensitys

def getTargets(mq, intensitys, mqs, maxDelQM, minIndensity):
	targets = [[filter(lambda elem: elem[0]>=minIndensity and elem[1]<mq00+maxDelQM and elem[1]>mq00-maxDelQM, zip(intensitys, mqs)) for mq00 in mq0] for mq0 in mq]
	targets = [filter(lambda elem: len(elem)>0, elem) for elem in targets]
	targets = [[reduce(lambda x,y: x if x>y else y, elem) for elem in target] for target in targets]
	return targets

def plotSpectra(bTargets, yTargets, mqs, intensitys, saveFigName):
	import matplotlib.pyplot as mpy
	fig, ax = mpy.subplots()
	markerline0, stemlines0, baseline0 = ax.stem(mqs, intensitys, markerfmt=' ')
        mpy.setp(stemlines0, color='k')
	for bTarget in bTargets:
		plotX, plotY = [], []
	        for y,x in bTarget:
			plotX.append(x)
	                plotY.append(y)
        	markerline1, stemlines1, baseline1 = ax.stem(plotX, plotY, markerfmt=' ')
	        mpy.setp(stemlines1, color='g')
	for yTarget in yTargets:
        	plotX, plotY = [], []
	        for y,x in yTarget:
        	        plotX.append(x)
                	plotY.append(y)
	        markerline2, stemlines2, baseline2 = ax.stem(plotX, plotY, markerfmt=' ')
        	mpy.setp(stemlines2, color='m')
	fig.show()
	fig.savefig(saveFigName)

def main():
	import sys
	massFileName = 'mass.ini'
	modifyFileName = 'modify.ini'
	maxDelQM = 0.5
	minIndensity = 0.01
	saveFigName = 'resultFig.eps'
	modifyIndexes=[]
	modifyNames=[]
	if len(sys.argv)>1:
		spectraFileName = sys.argv[1]
		peptide = sys.argv[2]
		charge = int(sys.argv[3])
	if len(sys.argv)>4:
		saveFigName = sys.argv[4]
	if len(sys.argv)>5:
		maxDelQM = float(sys.argv[5])
		minIndensity = float(sys.argv[6])
	if len(sys.argv)>7:
		modifyIndexes = [int(elem) for elem in sys.argv[7].split(',')]
		modifyNames = sys.argv[8].split(',')
	if len(sys.argv)>9:
		massFileName = sys.argv[9]
		modifyFileName = sys.argv[10]
	print 'massFileName='+massFileName+'\n'
	print 'modifyFileName='+modifyFileName+'\n'
        print 'spectraFileName='+spectraFileName+'\n'
        print 'peptide='+peptide+'\n'
        print 'charge='+str(charge)+'\n'
	if not len(modifyIndexes)==0:
	        print 'modifyIndexes='+reduce(lambda x,y: str(x)+','+str(y), modifyIndexes)+'\n'
	if not len(modifyNames)==0:
	        print 'modifyNames='+",".join(modifyNames)+'\n'
        print 'maxDelQM='+str(maxDelQM)+'\n'
        print 'minIndensity='+str(minIndensity)+'\n'
	print 'saveFigName='+saveFigName+'\n'
	massMap = getMassMap(massFileName)
	modifyMap = getModifyMap(modifyFileName)
	b, y = getPeptideMQ(peptide, charge, modifyIndexes, modifyNames, massMap, modifyMap)
	mqs, intensitys = getSpectraData(spectraFileName)
	bTargets = getTargets(b, intensitys, mqs, maxDelQM, minIndensity)
	yTargets = getTargets(y, intensitys, mqs, maxDelQM, minIndensity)
	plotSpectra(bTargets, yTargets, mqs, intensitys, saveFigName)

if __name__ == 	'__main__':
	main()
	
	
	





