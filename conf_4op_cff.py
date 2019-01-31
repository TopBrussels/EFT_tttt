config={'annotation': 'Example annotation',
 'mode': 'jinja-test',
 'combine_asymptotic_limits_json':'SL_OS_SS_combined_blind.json',
    'tables':{
        'independent':{
            'latex_main':'build/table_independent.tex',
            'template':'resources/common/table_independent_4op.jinja2.tex',
            'operators':['O_R','O_L^1','O_B^1','O_B^8']
        },
        'marginal':{
            'latex_main':'build/table_marginal.tex',
            'template':'resources/common/table_marginal_4op.jinja2.tex',
            'operators':['O_R','O_L^1','O_B^1','O_B^8']
        }
    }
 }
