#!/usr/bin/env python

"""
A simple python script template.
"""
__version__ = "1.0"

import os
import sys
import imp
import time
import glob
import shutil
import argparse
import subprocess
import logging
import json
import textwrap

import pprint as pp
import ROOT as rt

import jinja2
from jinja2 import Template

from eft_coefficients import EftPredictions
from mg_calculations import wilson_coefficients, MG_SM, sig_SM

import functools, logging

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class log_with(object):
    '''Logging decorator that allows you to log with a specific logger.'''
    # Customize these messages
    ENTRY_MESSAGE = 'Entering {}'
    EXIT_MESSAGE = 'Exiting {}'

    def __init__(self, logger=logging):
        self.logger = logger

    def __call__(self, func):
        '''Returns a wrapper that wraps func. The wrapper will log the entry and exit points of the function with logging.INFO level.'''
        # set logger if it was not set earlier
        if not self.logger:
            logging.basicConfig()
            self.logger = logging.getLogger(func.__module__)
            # self.logger.setLevel(logging.DEBUG)

        @functools.wraps(func)
        def wrapper(*args, **kwds):
            self.logger.debug(self.ENTRY_MESSAGE.format(func.__name__))  # logging level .info(). Set to .debug() if you want to
            f_result = func(*args, **kwds)
            self.logger.debug(self.EXIT_MESSAGE.format(func.__name__))   # logging level .info(). Set to .debug() if you want to
            return f_result
        return wrapper

def progress(current, total, status=''):
        fullBarLength = 80
        doneBarLength = int(round(fullBarLength * current / float(total)))

        percents = round(100.0 * current / float(total), 1)
        bar = '>' * doneBarLength + ' ' * (fullBarLength - doneBarLength)

        sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
        sys.stdout.flush()

class Model(object):
        def __init__(self,configuration):
                self._configuration = configuration
                self._sigma_hat = {}
                self._sigma_hat_ij = {}
                self._annotation = 'Performance comparision of different MVA discriminants'
                if 'annotation' in self._configuration:
                        self._annotation = self._configuration['annotation']
                self._eft = None
                self.Initialize()

        @log_with()
        def Initialize(self):
                self._eft = EftPredictions(wilson_coefficients, MG_SM, sig_SM)
                pass

        @log_with()
        def get_linear_coefficient(self,operator):
                """
                Factory method
                """
                if operator in self._sigma_hat:
                        return self._sigma_hat[operator]
                else:
                        if operator == 'O_R':
                                self._sigma_hat[operator] = self._eft.Sigma1_vec[0]
                        elif operator == 'O_L^1':
                                self._sigma_hat[operator] = self._eft.Sigma1_vec[1]
                        elif operator == 'O_L^8':
                                self._sigma_hat[operator] = self._eft.Sigma1_vec[2]
                        elif operator == 'O_B^1':
                                self._sigma_hat[operator] = self._eft.Sigma1_vec[3]
                        elif operator == 'O_B^8':
                                self._sigma_hat[operator] = self._eft.Sigma1_vec[4]
                        else: raise RuntimeError("Operator {} is not available".format(bcolors.BOLD+operator+bcolors.ENDC))
                return self._sigma_hat[operator]

        @log_with()
        def get_quadratic_coefficient(self, pair_of_operators):
                """
                Factory method
                """
                if len(pair_of_operators) != 2:
                        raise RuntimeError("Operator pair {} must have two entries".format(bcolors.BOLD+operator+bcolors.ENDC))
                op1, op1_ind     = pair_of_operators[0], self._eft.name_index_map.get(pair_of_operators[0],None)
                op2, op2_ind     = pair_of_operators[1], self._eft.name_index_map.get(pair_of_operators[1],None)
                if op1_ind is None or op2_ind is None:
                        raise RuntimeError("Operator {} or {} not found".format(bcolors.BOLD+op1+bcolors.ENDC, 
                                                                                bcolors.BOLD+op2+bcolors.ENDC))

                if op1 in self._sigma_hat_ij and op2 in self._sigma_hat_ij[op1]: return self._sigma_hat_ij[op1][op2]
                else: 
                        if op1 not in self._sigma_hat_ij: 
                                # print "1st case: ", op1,  self._sigma_hat_ij
                                self._sigma_hat_ij[op1]={op2:self._eft.Sigma2_matr[op1_ind][op2_ind]}
                        else: 
                                self._sigma_hat_ij[op1][op2]=self._eft.Sigma2_matr[op1_ind][op2_ind]

                return self._sigma_hat_ij[op1][op2]

class Serializer(object):
        @log_with()
        def __init__(self,builddir='build'):
                self._buildfolder = builddir
                self._outputfolder = None
                pass
        
        @log_with()
        def set_outputfolder(self,folder):
                self._outputfolder = folder
                if not os.path.exists(folder):
                        os.makedirs(folder)

        @log_with()
        def move_builddir_to_outputfolder(self):
                print self._buildfolder, self._outputfolder, (self._buildfolder and self._outputfolder)
                if self._buildfolder is not None and self._outputfolder is not None:
                        for extension in ['pdf','png','tex']:
                                for file in glob.glob('{}/*.{}'.format(self._buildfolder,extension)):
                                        shutil.move(file, self._outputfolder)

        @log_with()
        def serialize_view(self,View):
                self.move_builddir_to_outputfolder()
                pass
        
        @log_with()
        def serialize_beamer_view(self,View):
                self.move_builddir_to_outputfolder()
                pass

        
        @log_with()
        def serialize_report_view(self,View):
                self.move_builddir_to_outputfolder()
                pass

class View(object):
        @log_with()
        def __init__(self):
                self.model = None
                self._outfilename = 'out'
                self._outfileextension = 'png'
                self._outputfolder = 'build'
        @log_with()
        def set_model(self,model):
                self.model = model
        @log_with()
        def set_builddir(self,folder):
                self._outputfolder = folder
                if not os.path.exists(folder):
                        os.makedirs(folder)
        @log_with()
        def set_outfilename(self,filename):
                if filename: self._outfilename = filename
        @log_with()
        def set_extension(self,extension):
                self._outfileextension = extension
        @log_with()
        def get_outfile_name(self,substring=''):
                for ext in self._outfileextension.split(","):
                        yield '{}/{}{}.{}'.format(self._outputfolder,self._outfilename,substring,ext)

        @log_with()
        def annotate(self,type):
                if type == "screen":
                        bright_green_text = "\033[1;32;40m"
                        normal_text = "\033[0;37;40m"
                        print "\n".join(textwrap.wrap(bcolors.OKBLUE+
                                  self.model._annotation.encode('ascii')+
                                  bcolors.ENDC, 120))
                elif type == "tex":
                        logging.warning("Annotation format: {}. Not implemented yet!".format(type))
                elif type == "md":
                        logging.warning("Annotation format: {}. Not implemented yet!".format(type))
                else:
                        logging.error("Annotation format not recognized: {}".format(type))

        @log_with()
        def save_config(self, config):
                if os.path.exists(self._outputfolder):
                        # Writing configuration data
                        if "_cff.py" in config: 
				with open(self._outputfolder+'/'+os.path.basename(config), 'w') as f:
                                        serialized_config_str = pp.pformat(self.model._configuration)
                                        serialized_config_str = 'config='+serialized_config_str
                                        f.write(serialized_config_str)
                        elif ".json" in config: 
                                with open(self._outputfolder+'/'+os.path.basename(config), 'w') as f:
                                        json.dump(self.model._configuration, f, indent=4, sort_keys=True)

        @log_with()
        def save(self,serializer):
                serializer.serialize_view(self)

        @log_with()
        def draw(self):
                pass

class LatexBeamerView(View):
        @log_with()
        def __init__(self):
                self.limits=None
                pass

        @log_with()
        def Init(self):
                pass

        @log_with()
        def save(self,serializer):
                serializer.serialize_beamer_view(self)

        @log_with()
        def draw(self):
                self.Init()
		View.draw(self)
                print self.model._configuration
                subprocess.call(["pdflatex", "-interaction=nonstopmode", "-output-directory={}".format(self._outputfolder), 
                                 self.model._configuration['latex_main']])

class LatexReportView(View):
        @log_with()
        def __init__(self):
                self.latex_jinja_env = jinja2.Environment(
                        block_start_string = '\BLOCK{',
                        block_end_string = '}',
                        variable_start_string = '\VAR{',
                        variable_end_string = '}',
                        comment_start_string = '\#{',
                        comment_end_string = '}',
                        line_statement_prefix = '%%',
                        line_comment_prefix = '%#',
                        trim_blocks = True,
                        autoescape = False,
                        loader = jinja2.FileSystemLoader(os.path.abspath('.'))
                )

        @log_with()
        def Init(self):
                pass


        @log_with()
        def save(self,serializer):
                serializer.serialize_report_view(self)

        @log_with()
        def draw(self):
                self.Init()
		View.draw(self)
                print self.model._configuration
                tables_templates = {}
                for tab_name, tab_config in self.model._configuration['tables'].iteritems():
                        sigmas = {}
                        for entry in tab_config['operators']:
                                sigma_hat_i = None
                                sigma_hat_ij = None
                                if 'linear' in tab_name:
                                        sigma_hat_i = self.model.get_linear_coefficient(entry)
                                        sigmas[entry] = round(sigma_hat_i,2)
                                elif 'quadratic' in tab_name:
                                        # raise RuntimeError("Table {} is not supported. Terminating!".format(bcolors.BOLD+tab_name+bcolors.ENDC))
                                        sigma_hat_ij = self.model.get_quadratic_coefficient(entry)
                                        if entry[0] not in sigmas: sigmas[entry[0]] = {entry[1]:round(sigma_hat_ij,2)}
                                        else: sigmas[entry[0]][entry[1]] = round(sigma_hat_ij,2)
                                else: 
                                        raise RuntimeError("Table name {} is not recognized. Terminating!".format(bcolors.BOLD+tab_name+bcolors.ENDC))
                        
                        # Render and save table
                        self.tab_template = self.latex_jinja_env.get_template(tab_config['template'])
                        print sigmas
                        table_template = self.tab_template.render(sigma_hat=sigmas)
                        # determine output .tex file name by removing path and .jinja2 suffix
                        with open(tab_config['latex_main'], "w") as f:  # saves tex_code to outpout file
                                f.write(table_template)
                        tables_templates[tab_name] = table_template

                # Render and save report
                self.rep_template = self.latex_jinja_env.get_template(self.model._configuration['latex_main_template'])
                for tab_name, tab_config in self.model._configuration['tables'].iteritems():
                        logging.debug(tables_templates[tab_name])
                        report_template = self.rep_template.render(abstract=self.model._annotation,**tables_templates)
                        # determine output .tex file name by removing path and .jinja2 suffix
                        report_file_name = self.model._configuration['latex_main_template'].split('/')[-1].replace('.jinja2','')
                        with open("build/"+report_file_name, "w") as f:  # saves tex_code to outpout file
                                f.write(report_template)
                #build pdf file
                subprocess.call(["pdflatex", "-interaction=nonstopmode", "-output-directory={}".format(self._outputfolder),
                                 "build/"+report_file_name])
		
def main(arguments):

        # Disable garbage collection for this list of objects
        rt.TCanvas.__init__._creates = False
        rt.TFile.__init__._creates = False
	rt.TH1.__init__._creates = False
	rt.TH2.__init__._creates = False
        rt.THStack.__init__._creates = False
        rt.TGraph.__init__._creates = False
        rt.TMultiGraph.__init__._creates = False
        rt.TList.__init__._creates = False
        rt.TCollection.__init__._creates = False
        rt.TIter.__init__._creates = False

        #Load configuration .json file
        configuration = None
	if ".json" in arguments.config:
        	with open(arguments.config) as json_data:
                	configuration = json.load(json_data)
                	logging.debug(pp.pformat(configuration))
	elif "_cff.py" in arguments.config:
                configuration_module = imp.load_source('my_config', arguments.config)
		configuration = configuration_module.config
		logging.debug(pp.pformat(configuration))

        model = Model(configuration)

        view = None
        if configuration['mode'] == 'beamer':
                print "beamer option is not implemented!"
                view = LatexBeamerView()
        elif configuration['mode'] == 'report':
                view = LatexReportView()
        else:
                view = View()
        view.set_model(model)
        view.set_builddir(arguments.builddir)
        view.set_outfilename(arguments.outfile)
        view.set_extension(arguments.extension)
        view.draw()
        serializer = Serializer(builddir=arguments.builddir)
        serializer.set_outputfolder(arguments.dir)
        # view.save(serializer)
	
	configuration['command']=' '.join(sys.argv)
        if arguments.annotation_format:
                view.annotate(arguments.annotation_format)
                view.save_config(arguments.config)

        return 0


if __name__ == '__main__':
        start_time = time.time()

        parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
        parser.add_argument('-o', '--outfile', help="Output file", default='test')
        parser.add_argument('-e', '--extension', help="Plot file extension (.C, .root, .png, .pdf)", default='png')
        parser.add_argument('--builddir', help="Build directory", default='build')
        parser.add_argument('--dir', help="Result output directory", default='.')
        parser.add_argument('-c', '--config', help=".json or _cff.py configuration file", required=True)
        parser.add_argument('-a', '--annotation_format', default="screen",\
                            help="Print annotation in given format (screen, tex, md)")
        parser.add_argument('--no-annotation', dest='annotation_format', action='store_false',\
                                                help="Disable annotating")
        parser.add_argument('-b', help="ROOT batch mode", dest='isBatch', action='store_true')
        parser.add_argument(
                        '-d', '--debug',
                        help="Print lots of debugging statements",
                        action="store_const", dest="loglevel", const=logging.DEBUG,
                        default=logging.WARNING,
                        )
        parser.add_argument(
                        '-v', '--verbose',
                        help="Be verbose",
                        action="store_const", dest="loglevel", const=logging.INFO,
                        )

        args = parser.parse_args(sys.argv[1:])

        print(args)

        logging.basicConfig(level=args.loglevel)

        logging.info( time.asctime() )
        exitcode = main(args)
        logging.info( time.asctime() )
        logging.info( 'TOTAL TIME IN MINUTES:' + str((time.time() - start_time) / 60.0))
        sys.exit(exitcode)
