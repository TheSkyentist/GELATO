#! /usr/bin/env python

""" Turn parameter file into LaTeX table """

# Packages
import copy
import argparse

# gelato supporting files
import gelato.ConstructParams as CP
import gelato.AdditionalComponents as AC

# Turn parameters into TeX file
def TeX(params):

    # Add premble
    tex = '\documentclass{article}\n\\usepackage{multirow}\n\\usepackage{amsmath}\n\\begin{document}\n\\begin{center}\n\\noindent\n\\begin{tabular}{c|cc|ccccc}\n\hline\hline\hline\n'
    # Header
    tex += '\multicolumn{5}{c|}{\\textbf{Species (Tie \\boldmath$z$ \& $\sigma$)}} & \multicolumn{3}{c}{\\textbf{Groups}}\\\\\n'
    # Subheader
    tex += '\\textbf{Name} & \multicolumn{1}{|c}{\\textbf{Line [\\AA]}} & \multicolumn{1}{c|}{\\textbf{Ratio}} & \\textbf{+Comp} & \multicolumn{1}{c|}{\\textbf{$\\to$Group}} & \\textbf{Name} & \\textbf{Tie \\boldmath$z$?} & \\textbf{Tie \\boldmath$\sigma$?} \\\\\hline\hline\n'

    # Iterate over Groups
    for i,group in enumerate(params['EmissionGroups']):
    
        # Number of lines in the group
        Glines = 0
        for s in group['Species']:
            for l in s['Lines']: Glines += 1
        
        # If empty group, add it
        if Glines == 0: 

            tex += 5*' & '+' & '.join(['\multirow{1}{*}{'+x+'}' for x in [group['Name'],slashbool(group['TieRedshift']),slashbool(group['TieDispersion'])]])
            tex += '\\\\\n'
            if i < len(params['EmissionGroups']):
                tex += '\cline{0-6}\n'

        # Multirow
        for j,species in enumerate(group['Species']):
            
            # Number of lines in the species
            Slines = 0
            for l in species['Lines']: Slines += 1
            if Slines == 0: Slines = 1 # If no lines, still want some width to row

            # Iterate over line
            for k,line in enumerate(species['Lines']):

                # If first species, add it
                if k == 0:
                    
                    # Add in Species Name
                    tex += ' & '.join(['\multirow{'+str(Slines)+'}{*}{'+str(x)+'}' for x in [species['Name']]])

                    # Add in Wavelengths and Relstrenghs
                    tex += '&' + str(line['Wavelength']) + ' & ' + str(line['RelStrength']).replace('None','-') + '&'

                    # Add in Components
                    # Get components 
                    components = '{\\footnotesize '+', '.join([AC.ComponentName(i) for i,f in enumerate(bin(species['Flag'])[2:][::-1]) if f == '1'])+'}'
                    flaggroups =  '{\\footnotesize '+', '.join(species['FlagGroups'])+'}'
                    tex += ' & '.join(['\multirow{'+str(Slines)+'}{*}{'+str(x)+'}' for x in [components,flaggroups]]) +'&'

                    # If first group, add it
                    if j == 0:
                        tex += ' & '.join(['\multirow{'+str(Glines)+'}{*}{'+x+'}' for x in [group['Name'],slashbool(group['TieRedshift']),slashbool(group['TieDispersion'])]])
                
                else: 
                    tex += '&' + str(line['Wavelength']) + ' & ' + str(line['RelStrength']).replace('None','-')

                # Add new line
                tex += '\\\\\n'

                # Add horizontal lines
                if k == Slines - 1:
                    if j == len(group['Species'])-1:
                        if i < len(params['EmissionGroups']):
                            tex += '\cline{0-7}\n'
                    else: tex += '\cline{0-4}\n'
                else:
                    tex += '\cline{2-3}\n'

    tex += '\hline\hline\hline\n\end{tabular}\n\\end{center}\n\\end{document}'
    return tex

def slashbool(b):
    if b:
        return '\\textbf{True} $\\mid$ False'
    return 'True $\\mid$ \\textbf{False}'

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