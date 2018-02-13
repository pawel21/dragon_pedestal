import numpy as np


class Event:
    nch=8
    roisize = 40
    size4drs=4*1024

    def __init__(self, roisize):
        self.roisize = roisize
        self.header = 0
        self.dms = 0
        self.fcntpps = 0
        self.cnt10mhz = 0
        self.evtcnt = 0
        self.trigcnt = 0
        self.cnt133mhz = 0
        self.dml = 0
        self.firstcap = np.zeros((self.nch, 2), dtype=np.int)
        self.lasttime = np.zeros((self.nch, 2, self.size4drs), dtype=np.int)
        self.flags = []
        self.samples = np.zeros((8, 2, roisize))

    def read(self, file):
        self.reset_time()
        data = file.read()
        self.header =  int.from_bytes(data[:64], byteorder='big')
        self.dms = int.from_bytes(data[:2], byteorder='big')
        self.cntpps = int.from_bytes(data[2:4], byteorder='big')
        self.cnt10mhz = int.from_bytes(data[4:8], byteorder='big')
        self.evtcnt = int.from_bytes(data[8:12], byteorder='big')
        self.trigcnt = int.from_bytes(data[12:16], byteorder='big')
        self.cnt133mhz = int.from_bytes(data[16:24], byteorder='big')
        self.dml = int.from_bytes(data[24:32], byteorder='big')

        for i in range(0, 8):
            self.flags.append(data[32 + i * 2: 34 + i * 2])

        for i in range(0, self.nch):
            ich = int(i / 2) * 2
            gain = int(i % 2)
            value = int.from_bytes(data[48 + i * 2: 50 + i * 2], byteorder="big")
            self.firstcap[ich, gain] = value
            self.firstcap[ich + 1, gain] = value

        for k in range(0, 2 * self.roisize):
            d = data[(k * 16) + 64: 16 * (k + 1) + 64]
            for i in range(0, 4):
                if k < self.roisize:
                    pix = i * 2
                else:
                    pix = i * 2 + 1
                self.samples[pix, 0, k % self.roisize] = (int.from_bytes(d[i * 4: i * 4 + 2], byteorder="big"))
                self.samples[pix, 1, k % self.roisize] = (int.from_bytes(d[i * 4 + 2: i * 4 + 4], byteorder="big"))

        self.corr_time()

    def reset_time(self):
        for i in range(0, self.nch):
            for j in range(0, 2):
                for k in range(0, self.size4drs):
                    self.lasttime[i, j, k]=-1
                self.firstcap[i, j] = 0

    def get_DRS_time(self):
        return self.cnt133mhz*(1./133.e6)*1.e3;

    def corr_time(self):
        print("Corr time !!!!!!")
        timenow = self.get_DRS_time();
        for i in range(0, self.nch):
            for j in range(0, 2):
                fc = self.firstcap[i][j]
                #print(fc)
                for k in range(0, self.roisize):
                    posabs = (k+fc)%self.size4drs
                    if (self.lasttime[i, j, posabs]>0):
                        timediff = timenow - self.lasttime[i, j, posabs]
                        timecorr = self.pedtimecorr(timediff)
                        self.samples[i,j,k] -= timecorr
                    if k<self.roisize-1:
                        self.lasttime[i, j, posabs] = timenow
                    if i%2 == 0:
                        if fc%1024>766 and fc%1024 < 1012:
                            for kk in range(fc+1024-1, fc+1024+11):
                                self.lasttime[i, j, int(kk%self.size4drs)] = timenow
                        elif fc%1024 >= 1012:
                            channel = int(fc/1024)
                            for kk in range(fc+1024, (channel+2)*1024):
                                self.lasttime[i, j, int(kk%4096)] = timenow

    @staticmethod
    def pedtimecorr(timediff):
        if timediff > 0 and timediff<45:
            return 29.3*(timediff**(-0.2262))-12.4
        else:
            return 0

try:
    f = open("Randome7kHz20kev_run1.dat", "rb")
    ev = Event(40)
    ev.read(f)
except Exception as err:
    print(err)
finally:
    f.close()
