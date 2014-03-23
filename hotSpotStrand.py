import winsound

import pickle
#import matplotlib.pyplot as plt

class hotSpots():
    def __init__(self,binS = 5000000):
        try:
            self.d = pickle.load( open( "save.p", "rb" ) ) 
            print 'loaded'
        except:
            print 'reading GTF'
            self.d = self.makeGTFdict('Mus_musculus.GRCm38.74.gtf')
            pickle.dump(self.d, open( "save.p", "wb" ) )
        self.chrSize=dict.fromkeys(self.d.keys())
        for CHR in self.chrSize:
            for gene in self.d[CHR]:
                if self.chrSize[CHR] < self.d[CHR][gene][1]: self.chrSize[CHR]=self.d[CHR][gene][1]
        self.values = {'inGenes':1,'strandess':1}
        self.binS = binS
        self.b={}
        print 'bining for bins size: ' + str(binS)

    def getGTFline(self, line):
        details = line.split('\t')
        name=details[-1].split('gene_name "')[-1]
        name = name[:name.find('"')]
        chromosome =details[0]
        TYPE =details[1]
        dictionary={'chr':chromosome,'type':TYPE,'start':details[3],'end':details[4],'strand':details[6],'name':name}
        return dictionary
    
    def inGenesScore(self, CHR, start, end):
        Bin = 1
        d = self.d
        graph = dict.fromkeys(range(int(start/Bin), int(end/Bin)),0)
        genes = []
        for gene in d[CHR]:
            s = d[CHR][gene][0]/Bin
            e = d[CHR][gene][1]/Bin
            for i in range(s,e):
                try:
                    if not gene in genes: genes.append(gene)
                    graph[i]+=1
                    print 'gStart: ' + str(start/Bin) +' | gEnd: ' + str(end/Bin) + ' | index: ' + str(i) + ' | gene: ' +  gene
                except: pass
        score = 0
        prev=0
        
        keys = []
        counts = []
        for i in graph:
            keys.append(i)
            counts.append(graph[i])
            if graph[i]>1:
                if not prev == graph[i]:
                    prev = graph[i]
                    score+=graph[i]-1
            else: prev = 0
        

        #plt.plot(keys, counts)
        #plt.show()
        return score
        
    def genes(self, CHR, start, end):
        d = self.d
        tally = {'+':0,'-':0}        
        for gene in d[CHR]:
            a = d[CHR][gene][0]
            A = start
            b = d[CHR][gene][1]
            B = end
            x=0
            if a<A: x+=1
            if a<B: x+=1
            if b<A: x+=1
            if b<B: x+=1
            
            if x in [1,2,3]:
                tally[d[CHR][gene][2]]+=1
                
                #print 'start: ' +str(d[CHR][gene][0]-start)+' end: '+str(d[CHR][gene][1]-end)
        try:score = sum([tally['+'],tally['-']])
        except: score = 0
        return score
        
    def standessGenesScore(self, CHR, start, end):
        d = self.d
        tally = {'+':0,'-':0}

        
        for gene in d[CHR]:
            a = d[CHR][gene][0]
            A = start
            b = d[CHR][gene][1]
            B = end
            x=0
            if a<A: x+=1
            if a<B: x+=1
            if b<A: x+=1
            if b<B: x+=1
            
            if x in [1,2,3]:
                tally[d[CHR][gene][2]]+=1
                
                #print 'start: ' +str(d[CHR][gene][0]-start)+' end: '+str(d[CHR][gene][1]-end)
        try:score = 2 * min([tally['+'],tally['-']])/float(sum([tally['+'],tally['-']]))
        except: score = 0
        return score
    def makeGTFdict(self, FILE, allowedTypes = ['protein_coding']):
        f= open(FILE,'r')
        bothDicts={}
        names=[]
        while True:
            line = f.readline()
            if line == '':
                break
            line = self.getGTFline(line)
            if not(line['name'] in names):
                names.append(line['name'])
                if not line['chr'] in bothDicts: 
                    bothDicts.update({line['chr']:{}})
                    print line['chr']
                if line['name'] in bothDicts[line['chr']]:
                    if int(line['start']) < bothDicts[line['chr']][line['name']][0]:
                        bothDicts[line['chr']][line['name']][0]=int(line['start'])
                    if int(line['end']) > bothDicts[line['chr']][line['name']][1]:
                        bothDicts[line['chr']][line['name']][1]=int(line['end'])
                else:
                    bothDicts[line['chr']].update({line['name']:[int(line['start']),int(line['end']),line['strand']]})
            
        f.close()
        return bothDicts
    def bins(self):
        binSize = self.binS
        gftDict = self.d
        binsD = dict.fromkeys(gftDict.keys())
        for i in binsD: binsD[i]={}
        for CHR in binsD:

            for gene in gftDict[CHR]:
                start, end= gftDict[CHR][gene][0]/binSize, gftDict[CHR][gene][1]/binSize
                for i in range(start,end+1):
                    try: binsD[CHR][i*5000000]+=1
                    except: 
                        binsD[CHR].update({i*5000000:1})
        self.b = binsD
    def graphBins(self , show = False):
        binsD= self.b
        binS = self.binS
        plotD=dict.fromkeys(self.b.keys())
        for CHR in binsD:
            xAxis=sorted(binsD[CHR].keys())
            yAxis=[]
            for x in xAxis:
                yAxis.append(binsD[CHR][x])
            plotD[CHR]=(xAxis,yAxis)
            
        for CHR in plotD.copy():
            if len(plotD[CHR][0])<10:
                del plotD[CHR]
        n=0
        cols = 4
        f, axarr = plt.subplots(len(plotD)/cols+1,cols,figsize=(15, 15))
        for CHR in sorted(plotD.keys()):
            if len(plotD[CHR][0])>1:            
    
                axarr[n/cols,n%cols].scatter(plotD[CHR][0],plotD[CHR][1])
                axarr[n/cols,n%cols].set_title('chromosome ' + CHR)
    
                n+=1 
            else: pass
                
        plt.text(0.1,0.1, 'Bin size ' + str(binS))
        margin = 0.03
        top = 1-margin
        bottom = 0+margin
        left = 0+margin
        right = 1-margin
        hspace = 0.5
        wspace = 0.3
        plt.subplots_adjust(hspace=hspace,wspace=wspace, right = right, left = left, bottom = bottom, top = top)
        if show: plt.show()
        plt.savefig( str(binS)+'.png')
        return plotD
        
    def show(self):
        self.graphBins(True)
        
    def topRegions(self,n = 0):
        allEntries = []
        allEntriesKeys = []
        output = []
        for CHR in self.p:
            for i in range(len(self.p[CHR][0])):
                count = self.p[CHR][1][i]
                start = self.p[CHR][0][i]
                end = start+self.binS
                entry= [count,CHR,start,end]
                allEntriesKeys.append(count)
                allEntries.append(entry)
        
        top =  sorted(allEntriesKeys)[-n:][::-1]
        for i in allEntries:
            if i[0] in top:
                output.append(i)
        pickle.dump(output, open( "hotSpots"+str(self.binS)+".p", "wb" ) )
        self.writeDict(output)
        return top
        
            
    def writeDict(self, dictionary):
   	output = []
   	for i in dictionary:
  		output.append(str(i[0])+'\t'+i[1] + ':' + str(i[2]) + '-' + str(i[3]))
   	output = '\n'.join(output)
   	f = open('hotSpots'+str(self.binS) +'.csv','w')
   	f.write(output)
   	f.close()

    def windowBin(self, CHR, steps, window, strands=True):
        print 'windows binning'
        b={}
        start=0
        strandScore = 0
        inScore = 0
        d = self.d
        if not CHR in b:
            print
            length = self.chrSize[CHR]/steps
            b.update({CHR:dict.fromkeys(range(0,length,steps))})
        start=0
        while True:
            if strands: strandScore = self.standessGenesScore( CHR, start, start+window)
            genes = self.genes( CHR, start, start+window)
            score=[strandScore,genes]
            b[CHR][start]=score
            if start%(window*10) == 0:
                done = str(int(float(start)/self.chrSize[CHR]*100))
                print 'Chromosome : ' + CHR + ' | position: ' + str(start) + ' | done: ' + done + '% | score: ' + str(score)+' | inGene: '+str(inScore)+' | strandess: '+str(strandScore)
            if start > self.chrSize[CHR]: break
            start+=steps
        params = ''
        if strands: params = params+'_strands'
        name = 'chr%(CHR)sw%(window)s_s%(steps)s%(params)s.csv' % {'CHR':CHR,'window':window,'steps':steps,'params':params}
        self.b.update(b) 
        txt = []
        for i in b[CHR]:
            line = str(i)+'\t'+str(b[CHR][i][0])+'\t'+str(b[CHR][i][1])+'\n'
            txt.append(line)
        txt=''.join(txt)
        w = open(name,'w')
        w.write(txt)
        w.close() 
        print name

allCHR = ['JH584293.1',
 'GL456354.1',
 'MG4180_PATCH',
 'GL456219.1',
 'GL456381.1',
 'MG132_PATCH',
 'MG4136_PATCH',
 'GL456385.1',
 'JH584298.1',
 'MG4213_PATCH',
 'MG153_PATCH',
 'GL456372.1',
 'MG4151_PATCH',
 '3',
 '2',
 '5',
 '4',
 '7',
 'MG4211_PATCH',
 '9',
 'MG4222_MG3908_PATCH',
 'GL456211.1',
 'JH584296.1',
 'MG4209_PATCH',
 'MG3829_PATCH',
 'GL456221.1',
 'GL456350.1',
 '6',
 'GL456233.1',
 'MG4212_PATCH',
 '17',
 'MG4214_PATCH',
 '8',
 'Y',
 'X',
 'GL456239.1',
 'JH584294.1',
 '11',
 '10',
 '13',
 '12',
 '15',
 '14',
 'GL456210.1',
 '16',
 '19',
 '18',
 'JH584292.1',
 'JH584295.1',
 'GL456216.1',
 'MT',
 'JH584297.1',
 '1',
 'GL456212.1',
 'JH584304.1',
 'MG3833_PATCH',
 'JH584299.1',
 'JH584303.1']
 
hot = hotSpots()
#print hot.standessGenesScore('4',0,5000000)
for c in allCHR:
    hot.windowBin(c,50000,1000000)

#hot.show()
'''
hot.topRegions(100)
'''
winsound.Beep(1000,1000)
