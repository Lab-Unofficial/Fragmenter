

class Fragmenter():
    """Fragmenter takes a ssDNA sequence, generates a reverse complement, and appends restriction enzyme sequences to the 5' and 3' ends. This utility assumes that all cut site sequences are contained within the backbone segment being excised. Replace N nucleotides with the corresponding sequence in the plasmid backbone receiving the fragment.
    """

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    RE = {
        'NONE' : {'5P' : '',
                '3P' : '',
                'DIGEST': 'NULL'},

        'BBSI' : {'5P' : 'NNNN',
                '3P' : '',
                'DIGEST' : 'GAAGAC'},

        'SAPI' : {'5P' : 'NNN',
                '3P' : '',
                'DIGEST': 'GCTCTT'}
    }

    def __init__(self):
        self._seq = ''
        self._rev_comp = ''
        self._5_RE = 'NONE'
        self._3_RE = 'NONE'
        self._HH_5_RE = 'NONE'
        self._HH_3_RE = 'NONE'

    @property
    def sequence(self):
        return self.RE[self._5_RE]['5P']+ self._seq + self.RE[self._3_RE]['3P']

    @sequence.setter
    def sequence(self, seq: str):
        """Check that bases in seq are ATGC, store if so."""

        temp_seq = seq.upper()
        num_unknowns = sum([base not in self.complement for base in temp_seq])
        if num_unknowns:
            raise ValueError('A base was not recognized! Please use ATGC only. Sequence provided: {}'.format(seq))
        self._seq = temp_seq
        self._rev_comp = self._revComp()

    @property
    def reverse_complement(self):
        return self.RE[self._3_RE]['5P'] + self._rev_comp + self.RE[self._5_RE]['3P']

    @property
    def report(self):

        cut_conflict = [ self.RE[key]['DIGEST']  in self._seq + self._rev_comp  for key in self.RE.keys()]

        danger_zyme = []
        for idx, val in enumerate(cut_conflict):
            if val == True:
                danger_zyme.append(list(self.RE.keys())[idx])

        spacer_5 = (' ' * (len(self.RE[self.RE_5]['3P']) - len(self.RE[self.RE_5]['5P'])) )
        spacer_3 = (' ' * (len(self.RE[self.RE_5]['5P']) - len(self.RE[self.RE_5]['3P'])))

        HH_seqs = self._hammerhead()

        spacer_5_HH = (' ' * (len(self.RE[self.RE_5_HH]['3P']) - len(self.RE[self.RE_5_HH]['5P'])) )
        spacer_3_HH = (' ' * (len(self.RE[self.RE_5_HH]['5P']) - len(self.RE[self.RE_5_HH]['3P'])))

        return print('\nSequence provided: {}\n' \
                    'Fragment 5\u0027 RE: {}\n' \
                    'Fragment 3\u0027 RE: {}\n' \
                    'Forward fragment ({} nucleotides)\n' \
                    'Reverse complement fragment ({} nucleotides)\n\n' \

                    'fwd      5\u0027\u2794 3\u0027 : {}{}\n' \
                    'rev comp 3\u0027\u2794 5\u0027 : {}{}\n' \
                    'rev comp 5\u0027\u2794 3\u0027 : {}\n\n' \

                    'Hammerhead 5\u0027 RE: {}\n' \
                    'Hammerhead 3\u0027 RE: {}\n' \
                    'FWD HH fragment ({} nucleotides)\n' \
                    'REV HH fragment ({} nucleotides)\n\n' \

                    'HH fwd      5\u0027\u2794 3\u0027 : {}{}\n' \
                    'HH rev comp 3\u0027\u2794 5\u0027 : {}{}\n' \
                    'HH rev comp 5\u0027\u2794 3\u0027 : {}\n\n' \

                    'Cut sites found in fragments: {}\n' \
                    'NOTE: Replace all N nucleotides to match the plasmid backbone digest sites!\n'.format(self._seq, 
                                    self.RE_5, 
                                    self.RE_3, 
                                    len(self.sequence),
                                    len(self.reverse_complement),
                                    spacer_5, 
                                    self.sequence, 
                                    spacer_3, 
                                    self.reverse_complement[::-1],
                                    self.reverse_complement,
                                    self.RE_5_HH,
                                    self.RE_3_HH,
                                    len(HH_seqs[0]),
                                    len(HH_seqs[1]),
                                    spacer_5_HH,
                                    HH_seqs[0],
                                    spacer_3_HH,
                                    HH_seqs[1][::-1],
                                    HH_seqs[1],
                                    danger_zyme))

    @property
    def RE_5(self):
        return self._5_RE

    @RE_5.setter
    def RE_5(self, restriction_enzyme: str):
        temp_re = restriction_enzyme.upper()
        if temp_re in self.RE:
            self._5_RE = temp_re

        else:
            raise ValueError('Restriction enzyme not recognized! \nRE provided: {} \nREs accepted: {}'.format(restriction_enzyme, list(RE.keys())))

    @property
    def RE_5_HH(self):
        return self._HH_5_RE

    @RE_5_HH.setter
    def RE_5_HH(self, restriction_enzyme: str):
        temp_re = restriction_enzyme.upper()
        if temp_re in self.RE:
            self._HH_5_RE = temp_re

        else:
            raise ValueError('Restriction enzyme not recognized! \nRE provided:     {} \nREs accepted: {}'.format(restriction_enzyme, list(RE.keys())))

    @property
    def RE_3(self):
        return self._3_RE

    @RE_3.setter
    def RE_3(self, restriction_enzyme: str):
        temp_re = restriction_enzyme.upper()
        if temp_re in self.RE:
            self._3_RE = temp_re

        else:
            raise ValueError('Restriction enzyme not recognized! \nRE provided: {} \nREs accepted: {}'.format(restriction_enzyme, list(RE.keys())))

    @property
    def RE_3_HH(self):
        return self._HH_3_RE

    @RE_3_HH.setter
    def RE_3_HH(self, restriction_enzyme: str):
        temp_re = restriction_enzyme.upper()
        if temp_re in self.RE:
            self._HH_3_RE = temp_re

        else:
            raise ValueError('Restriction enzyme not recognized! \nRE provided: {} \nREs accepted: {}'.format(restriction_enzyme, list(RE.keys())))

    def _revComp(self):
        return ''.join([self.complement[base] for base in self._seq[::-1]])
        
    def _hammerhead(self):
        """Returns hammerhead complement of the last 10 nucleotides of the sequence. All returned elements are in the 5 to 3 direction. First element is for the forward strand and the second is for the reverse."""
        if len(self._seq) < 10:
            raise ValueError('Cannot generate hammerhead sequences. Fragment sequence must be >9 nucleotides.')
        if (self._HH_3_RE is not 'NONE') & (self._HH_5_RE is not 'NONE'):
            hh_fwd = self.RE[self._HH_5_RE]['5P'] + self._rev_comp[0:9] + self.RE[self._HH_3_RE]['3P']
            hh_rev_comp = self.RE[self._HH_3_RE]['5P'] + self._seq[-10:-1][::-1] + self.RE[self._HH_5_RE]['3P'] 
            return [ hh_fwd, hh_rev_comp]
        else:
            return ['','']
    
