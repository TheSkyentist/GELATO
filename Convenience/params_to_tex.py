#! /usr/bin/env python

""" Turn parameter file into LaTeX table """

# Packages
import copy
import argparse

# gelato supporting files
import gelato.ConstructParams as CP

# Turn parameters into TeX file
def TeX(params):

    tex = '\documentclass{article}\n\\usepackage{multirow}\n\\begin{document}\n\\begin{center}\n\\begin{tabular}{cccccccc}\n\hline\hline\hline\n'
    tex += '\\multicolumn{3}{c}{\\textbf{Groups}} & \multicolumn{3}{|c|}{\\textbf{Species}} & \multicolumn{2}{c}{\\textbf{Lines}}\\\\\n\\textbf{Name} & \\textbf{Tie \\boldmath$z$?} & \\textbf{Tie \\boldmath$\sigma$?} & \multicolumn{1}{|c}{\\textbf{Name}} & \\textbf{Flag} & \multicolumn{1}{c|}{\\textbf{Extra}} & \\textbf{\\boldmath$\lambda$ [\AA]} & \\textbf{Ratio}\\\\\hline\hline\n'

    for group in params['EmissionGroups']:
        # Number of lines
        Glines = 0
        for s in group['Species']:
            for l in s['Lines']: Glines += 1
        if Glines == 0: Glines = 1
        
        # Multirow
        tex += '&'.join(['\multirow{'+str(Glines)+'}{*}{'+x+'}' for x in [group['Name'],str(group['TieRedshift']),str(group['TieDispersion'])]]) + '&'
        
        for i,species in enumerate(group['Species']):
            Slines = 0
            for l in species['Lines']: Slines += 1
            if Slines == 0: Slines = 1

            if i > 0: tex += ' & & & '
            tex += ' & '.join(['\multirow{'+str(Slines)+'}{*}{'+str(x)+'}' for x in [species['Name'],species['Flag'],str(species['FlagGroups']).replace("'",'').replace('[','').replace(']','')]]) + '&'
            
            for j,line in enumerate(species['Lines']):

                if j > 0: tex += ' & & & & & &'
                tex += '' + str(line['Wavelength']) + ' & ' + str(line['RelStrength']).replace('None','-') + '\\\\\n'

            tex += '\cline{4-8}\n'

        if len(group['Species']) == 0: tex += '\\\\'
        tex += '\hline\hline\n'

    tex += '\hline\n\end{tabular}\n\\end{center}\n\\end{document}'
    return tex

# Main Function
if __name__ == "__main__":

    ## Parse Arguements to find Parameter File ##
    parser = argparse.ArgumentParser()
    parser.add_argument('Parameters', type=str, help='Path to parameters file')
    args = parser.parse_args()
    p = CP.construct(args.Parameters)
    ## Parse Arguements to find Parameter File ##

    ## Create Directory for Output
    if '.json' in args.Parameters:
        outname = args.Parameters.replace('.json','.tex')
    else: 
        outname = args.Parameters + '.tex'

    with open(outname,'w') as f:
        f.write(TeX(p))