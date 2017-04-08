"""
Contains a convenience dictionary relating the code quantities names
with their LaTeX representation as well as units.

Structure:
    quantities_dic = { code-name : ( LaTeX-name, LaTeX-units )}

Note: LaTeX-name and LaTeX units assume they are in the math mode!
"""

quantities_dict = {
  'h' : ( r'h', r'{\rm pc}'),
  'n' : ( r'n', r'{\rm cm}^{-3}'),
  'Bp' : ( r'\overline{ B }_\phi', r'\mu{\rm G}'),
  'Br' : ( r'\overline{ B }_r', r'\mu{\rm G}'),
  'Bzmod' : ( r'|\overline{ B }_z|', r'\mu{\rm G}'),
  'Beq' : ( r'\overline{ B }_{\rm eq}', r'\mu{\rm G}'),
  'Btot' : ( r'|\overline{\mathbf{ B } }|', r'\mu{\rm G}'),
  'Omega' : ( r'\Omega', r'{\rm km}\,{\rm s}^{-1}{\rm kpc}^{-1}'),
  'Shear' : ( r'S', r'{\rm km}\,{\rm s}^{-1}{\rm kpc}^{-1}'),
  'Uz' : ( r'U_z', r'{\rm km}\,{\rm s}^{-1}'),
  'p' : ( r'p', r'{\rm deg}'),
  'growth' : (r'\gamma', r'\rm Gyr^{{-1}}'),
  'Mstars_disk' : (r'M_{\star,\rm disk}', r'\rm M_\odot'),
  'Mstars_bulge' : (r'M_{\star,\rm bulge}', r'\rm M_\odot'),
  'Mgas_disk' : (r'M_{\rm g, disk}', r'\rm M_\odot'),
  'Bmax' : (r'\max(B)', r'\mu{\rm G}')
  }

