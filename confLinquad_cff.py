config={'annotation': 'Linear and quadratic dependeses of tttt cross section on the 4TopEFT degrees of freedom.',
 'command': 'linear_quadratic_coef.py -b -c confLinquad_cff.py',
 'mode': 'report',
 'latex_main_template': 'resources/report.jinja2.tex',
 'tables':{
    'table_linear':{
        'latex_main':'build/table_linear.tex',
        'template':'resources/common/table_tttt_eft_linear_dependence.jinja2.tex',
        'operators':['O_R','O_L^1','O_L^8','O_B^1','O_B^8']
    },
    'table_quadratic':{
        'latex_main':'build/table_quadratic.tex',
        'template':'resources/common/table_tttt_eft_quadratic_dependence_superdiagonal.jinja2.tex',
        'operators':[['O_R','O_R'],   ['O_R','O_L^1'],   ['O_R','O_L^8'],  ['O_R','O_B^1'],  ['O_R','O_B^8'],
                     ['O_L^1','O_R'], ['O_L^1','O_L^1'], ['O_L^1','O_L^8'],['O_L^1','O_B^1'],['O_L^1','O_B^8'],
                     ['O_L^8','O_R'], ['O_L^8','O_L^1'], ['O_L^8','O_L^8'],['O_L^8','O_B^1'],['O_L^8','O_B^8'],
                     ['O_B^1','O_R'], ['O_B^1','O_L^1'], ['O_B^1','O_L^8'],['O_B^1','O_B^1'],['O_B^1','O_B^8'],
                     ['O_B^8','O_R'], ['O_B^8','O_L^1'], ['O_B^8','O_L^8'],['O_B^8','O_B^1'],['O_B^8','O_B^8'],                     
        ]
    }
 }
}