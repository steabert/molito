                for line in open('elements.dat'):
                        if line.startswith('#'): continue
                        data = line.split()
                        self.PSE[data[1]] = Element(data[1], int(data[0]), float(data[3]), int(data[6]), numpy.array(data[11:14], dtype='float32'))
                        #my_file.write('\'{symbol}\':{radius}, '.format(symbol=data[1], radius=data[3]))
                        my_file.write(str(self.PSE[data[1]])+'\n')
