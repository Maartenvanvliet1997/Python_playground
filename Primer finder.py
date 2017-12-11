"""
seq = ['gacatactgagaataaatccaaagacattagtttctttgcacgaaatgaggttacatatccagtgacatttatttgagctatttaaacaacttaaacatctttttcttttcttaataagggacgtttcaagttgtggtctcagccaaa'
       ,'tttcaatatcagcacactcattctttgtcaattcattttttcccatgagatgaagcacatgtgacgaatacggactagataacctctaagaattttccacttcttcaaaatgaacttactctagaaagcttacccttggataaccagtttgactttcataatgtctctgttttttgtttttccaacaattacagactcaggttctcttattttggaagtttctatctggttttgttctgaacttacattttttttttttttggtatctatgattttttttgctcagggcatcaaaatgtgctaaggacaagaattatatcctttttaaaaaatgttgttagcttggtgtaaaatgtatattgactgtattggtgaataaattgaatagacataacctcaaagtacttcacttattctttttaactactgatttgataaaaagtatgattataagatatccacgacaatctcatagtttctt']
"""

p_nuc = {'A':'T','T':'A','G':'C','C':'G'}

class Primer:
    def __init__(self, name, Tm, GC, length):
        self.name = name
        self.Tm = Tm
        self.GC = GC
        self.length = length

    def get_name(self):
        return self.name

    def get_Tm(self):
        return self.Tm

    def get_GC(self):
        return self.GC

    def get_length(self):
        return self.length

class Application:
    def main(self):
        """
        a = input("Input 3' - 5' kant: ")
        b = input("Input 3' - 5' kant: ")
        seq = [a,b]
        """

        seq = ['gacatactgagaataaatccaaagacattagtttctttgcacgaaatgaggttacatatccagtgacatttatttgagctatttaaacaacttaaacatctttttcttttcttaataagggacgtttcaagttgtggtctcagccaaa',
               'tttcaatatcagcacactcattctttgtcaattcattttttcccatgagatgaagcacatgtgacgaatacggactagataacctctaagaattttccacttcttcaaaatgaacttactctagaaagcttacccttggataaccagtttgactttcataatgtctctgttttttgtttttccaacaattacagactcaggttctcttattttggaagtttctatctggttttgttctgaacttacattttttttttttttggtatctatgattttttttgctcagggcatcaaaatgtgctaaggacaagaattatatcctttttaaaaaatgttgttagcttggtgtaaaatgtatattgactgtattggtgaataaattgaatagacataacctcaaagtacttcacttattctttttaactactgatttgataaaaagtatgattataagatatccacgacaatctcatagtttctt']

        self.read_seq(seq)
        self.choose_primers()

    def read_seq(self, seq):
        self.f_primers = {}
        self.r_primers = {}

        for time in range(len(seq)):
            for i in range(len(seq[time])):
                primer = ''
                for char in seq[time][i:].upper():
                    if time:
                        primer += char
                    else:
                        primer += p_nuc[char]

                    self.is_primer(primer, time)

    def is_primer(self, primer, time):
        nuc = {'A':0, 'T':0, 'G':0, 'C':0}
        total = 0

        for char in primer:
            nuc[char] += 1
            total += 1

        Tm = (((nuc['A'] + nuc['T'])*2)+((nuc['G'] + nuc['C'])*4))
        GC = float(((nuc['G'] + nuc['C'])/total)*100)

        if primer[len(primer) - 1] not in 'GC':
            return

        if Tm < 55 or Tm > 62:
            return

        if GC < 40 or GC > 60:
            return

        if time:
            self.f_primers[primer] = Primer(primer, Tm, GC, total)
        else:
            self.r_primers[primer] = Primer(primer, Tm, GC, total)

    def choose_primers(self):
        all_primers = [self.f_primers, self.r_primers]
        with open('info_primers.txt', 'w') as f:
            for time in range(len(all_primers)):
                if time:
                    f.write('REVERSE PRIMERS:\n')
                else:
                    f.write('FORWARD PRIMERS:\n')
                for primer in all_primers[time]:
                    f.write(' - ' + primer +
                            ' Tm: ' + str(all_primers[time][primer].get_Tm()) +
                            ' GC: ' + ("%.2f" % all_primers[time][primer].get_GC()) + '%' +
                            ' Length: ' + str(all_primers[time][primer].get_length()) + '\n')


app = Application()
app.main()