from tkinter import *
import os

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
    def __init__(self):
        self.create_widgets()

    def create_widgets(self):
        #Creating widgets
        w = Tk()
        w.title('Primer finder')
        w.geometry('400x450')
        w.resizable(0, 0)

        l1 = Label(w, text="3' - 5' sequence: ")
        l2 = Label(w, text="5' - 3' sequence: ")

        self.u1 = StringVar()
        self.u2 = StringVar()

        e1 = Entry(w, textvariable=self.u1, width=35)
        e2 = Entry(w, textvariable=self.u2, width=35)

        b1 = Button(w, text='Process', command=self.start)
        b2 = Button(w, text='Show All', command=self.show_primers)

        l3 = Label(w, fg='blue', text="Forward Primer:")
        l4 = Label(w, fg='red', text="Reverse Primer:")

        self.f_p = StringVar()
        self.r_p = StringVar()

        self.f_p.set('X')
        self.r_p.set('X')

        l5 = Label(w, fg='blue', textvariable=self.f_p)
        l6 = Label(w, fg='red', textvariable=self.r_p)

        # Placing widgets on grid

        l1.grid(row=0, columnspan=2, sticky=W)
        l2.grid(row=1, columnspan=2, sticky=W)
        e1.grid(row=0, column=1)
        e2.grid(row=1, column=1)
        b1.grid(row=2, column=0)
        b2.grid(row=2, column=1)
        l3.grid(row=3, column=0)
        l4.grid(row=3, column=1)
        l5.grid(row=4, column=0)
        l6.grid(row=4, column=1)

        w.mainloop()
        ###############

    def start(self):
        self.read_seq()
        self.choose_primers()


    def read_seq(self):
        seq = [self.u1.get(), self.u2.get()[::-1]]

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

        chosen_primers = {}

        ###Calculates ^Tm

        for f_primer in self.f_primers:
            Tm_1 = self.f_primers[f_primer].get_Tm()
            for r_primer in self.r_primers:
                Tm_2 = self.r_primers[r_primer].get_Tm()

                Tm_diff = (max(Tm_1, Tm_2) - min(Tm_1, Tm_2))

                if Tm_diff <= 5:
                    if Tm_diff not in chosen_primers:
                        chosen_primers[Tm_diff] = [self.f_primers[f_primer], self.r_primers[r_primer]]
                    else:
                        len_1 = self.f_primers[f_primer].get_length()
                        len_2 = self.r_primers[r_primer].get_length()

                        if len_1 > chosen_primers[Tm_diff][0].get_length():
                            chosen_primers[Tm_diff][0] = self.f_primers[f_primer]

                        if len_2 > chosen_primers[Tm_diff][1].get_length():
                            chosen_primers[Tm_diff][1] = self.r_primers[r_primer]

        for obj in chosen_primers:
            self.f_p.set(chosen_primers[obj][0].get_name())
            self.r_p.set(chosen_primers[obj][1].get_name())
            return

    def show_primers(self):
        all_primers = [self.f_primers, self.r_primers]

        with open('info_primers.txt', 'w') as f:
            for time in range(len(all_primers)):
                if time:
                    f.write('*'*10 + 'REVERSE PRIMERS:' + '*' * 10 + '\n')
                else:
                    f.write('*' * 10 + 'FORWARD PRIMERS:' + '*' * 10 + '\n')
                for primer in all_primers[time]:
                    f.write(' - ' + primer + '\n' +
                            'Tm: ' + str(all_primers[time][primer].get_Tm()) +
                            ' GC: ' + ("%.2f" % all_primers[time][primer].get_GC()) + '%' +
                            ' Length: ' + str(all_primers[time][primer].get_length()) + '\n')

        os.startfile('info_primers.txt')

app = Application()

