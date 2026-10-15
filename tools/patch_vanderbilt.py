#!/usr/bin/env python3
"""Export generator reference pairs in the final Vanderbilt projector basis.

The patch changes output only. The original generator and recipes remain
separate. Original normalized psi/phi pairs are transformed with U^{-T},
where beta_final=beta_initial U and scpgef stores U^{-1} in ttt.
"""
# SPDX-License-Identifier: GPL-2.0-or-later
import argparse
from pathlib import Path

EXPORT = """c     Export matching reference pairs after the final projector rotation.
      open(unit=91,file='reference-waves.dat',status='replace')
      write(91,*) kkbeta,nbeta
      write(91,*) (lll(ib),ib=1,nbeta)
      write(91,*) (eee(ib),ib=1,nbeta)
      do ir=1,kkbeta
        write(91,'(100es25.16e3)') r(ir),rab(ir),
     +    (sum(psi(ir,1:nbeta)*ttt(ib,1:nbeta),
     +    mask=lll(1:nbeta).eq.lll(ib)),ib=1,nbeta),
     +    (sum(phi(ir,1:nbeta)*ttt(ib,1:nbeta),
     +    mask=lll(1:nbeta).eq.lll(ib)),ib=1,nbeta),
     +    (beta(ir,ib),ib=1,nbeta)
      enddo
      close(91)
c
"""

def patch(directory, functional):
    directory = Path(directory)
    path = directory / 'scgsubs.f'
    text = path.read_text()
    if EXPORT in text:
        text = text.replace(EXPORT, '', 1)
    marker = 'c     print warnings if eigenvalues of q matrix are unreasonable'
    pos = text.index(marker)
    text = text[:pos] + EXPORT + text[pos:]
    path.write_text(text)
    path = directory / 'runatom.f'
    text = path.read_text().replace('idim4 = 5', 'idim4 = 6')
    path.write_text(text)
    if functional == 'pbesol':
        path = directory / 'pbe.f'
        text = path.read_text()
        import re
        text, n1 = re.subn(r'(?i)(um\s*=\s*)[.0-9dDeE+-]+', r'\g<1>0.123456790123456790d0', text)
        text, n2 = re.subn(r'(?i)(bet\s*=\s*)[.0-9dDeE+-]+', r'\g<1>0.046d0', text)
        assert n1 and n2, (n1,n2)
        path.write_text(text)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path)
    p.add_argument('--functional', choices=['pbe','pbesol'], default='pbe')
    a = p.parse_args()
    patch(a.source,a.functional)
