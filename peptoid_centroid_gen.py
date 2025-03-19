# Peptoid.params to Peptoid_centroid.params script, 
# Allon Goldberg, Research Assistant, Flatiron Institute, 6/2024

import os
import re
import glob


def main(fa_inputs, output_dir):
    #Make dir if none exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    #Gather input files into vector according to given pattern
    fa_inputs = glob.glob(fa_inputs)
    #List of acceptable atoms in the params file for a peptoid centroid
    included_atoms = ['N', 'CA', 'C', 'O', 'CA1', 'LOWER', 'UPPER']
    for fa_input_file in fa_inputs:
        charge_line = ''
        #Set output file name as original name in new path
        filename = os.path.basename(fa_input_file)
        output_centroid_file = output_dir+'/'+filename
        #Raise error to prevent overwriting original files
        if os.path.abspath(fa_input_file) == os.path.abspath(output_centroid_file):
            raise ValueError('Will not overwrite original files, change output directory to new directory')

        #Read file
        with open(fa_input_file, 'r') as fa_params:
            lines = fa_params.readlines()
    
        header_index = 0
        line_number = 0
        cen_lines = []
        for line in lines:
            line_number += 1
            line_cache = line
            # Find where header ends
            if line_number <= 10 and line.startswith('#'):
                header_index += 1
            # Skip comments
            if line.startswith('#'):
                cen_lines.append(line)
                continue
            # Split line by arbitrary spaces (whitespace) into columns
            columns = re.split(r'\s+', line.strip())
            #Save name of peptoid
            if columns[0] == 'NAME':
                name = columns[1]
                cen_lines.append(line)
                continue
            #ATOM lines to be skipped
            if columns[0] == 'ATOM' and columns[1] not in included_atoms:
                continue
            elif columns[0] == 'ATOM' and columns[1] in included_atoms:
                cen_lines.append(line_cache)
                continue
            #BOND lines to be skipped
            if (columns[0] == 'BOND' or columns[0] == 'BOND_TYPE') and (columns[1] not in included_atoms or columns[2] not in included_atoms):
                continue
            #CHI lines to be skipped
            if columns[0] == 'CHI' or columns[0] == 'PROTON_CHI':
                continue
            #CHARGE line
            if columns[0] == 'CHARGE':
                columns[1] = 'CEN'
                charge_line = ' '.join(columns) + '\n'
                continue

            #NBR_ATOM line
            if columns[0] == 'NBR_ATOM':
                cen_lines.append('NBR_ATOM  CEN\n')
                continue
            #NBR_RADIUS line with radius adjustment
            if columns[0] == 'NBR_RADIUS':
                try:
                    columns[1] = str(round(float(columns[1])*0.536 + 2.33, 4))
                except ValueError:
                    continue
            #Store ICOOR line for centroid coordinates later
            if columns[0] == 'ICOOR_INTERNAL' and columns[1] == 'CA1':
                icoor_copy = line
            #ICOOR lines to be skipped or included
            if columns[0] == 'ICOOR_INTERNAL' and columns[1] not in included_atoms:
                continue
            elif columns[0] == 'ICOOR_INTERNAL' and columns[1] in included_atoms:
                cen_lines.append(line_cache)
                continue
            
            # Join columns back into a line with adjusted columns
            cen = ' '.join(columns) + '\n'
            cen_lines.append(cen)
        #Add centroid lines    
        cen_lines.append('##centroid-specific\n')
        if name == '207':
            AAeq = 'VAL'
        elif name == '501' or name == '502' or name == '504':
            AAeq = 'PHE'
        elif name == '004' or name == '006' or name == '010' or name == '111' or name == '018' or name == '019' or name == '021' or name == '108' or name == '110' or name == '113' or name == '132' or name == '404' or name == '405' or name == '407' or name == '413':
            AAeq = 'TYR'
        elif name == '503' or name == '505' or name == '506' or name == '127' or name == '128' or name == '129' or name == '130' or name == '009' or name == '020' or name == '005' or name == '411':
            AAeq = 'TRP'
        elif name == '002' or name == '008' or name == '011' or name == '012' or name == '014' or name == '015' or name == '112' or name == '131' or name == '020' or name == '005' or name == '411':
            AAeq = 'TRP'
        elif name == '507' or name == '508' or name == '509':
            AAeq = 'HIS'
        elif name == '701' or name == '702':
            AAeq = 'VAL'
        elif name == '703' or name == '704':
            AAeq = 'ILE'
        elif name == '313' or name == '314':
            AAeq = 'ASP'
        elif name == '332' or name == '333':
            AAeq = 'LYS'
        elif name == '601' or name == '623':
            AAeq = 'PHE'
        elif name == '621':
            AAeq = 'TRP'
        elif name[0] == '0' or name[0] == '1' or name[0] == '4':
            AAeq = 'PHE'
        elif name[0] == '2':
            AAeq = 'VAL'
        elif name[0] == '3' or name == '631' or name == '633':
            AAeq = 'SER'
        elif name == 'SAR':
            AAeq = 'ALA'
        else:
            AAeq = 'ALA'
        cen_lines.append('ATOM  CEN CEN_'+AAeq+' H 0.0\n')
        cen_lines.append('BOND N CEN\n')
        cen_lines.append(charge_line)
        #Maintain original spacing for centroid 'ball' coordinates with calculation
        icoor_split = re.split(r'(\s+)', icoor_copy.strip())
        icoor_split[8] = str(round(float(icoor_split[8])*0.544+0.483, 4))
        icoor_split[2] = 'CEN'
        icoor_cen_line = ''.join(icoor_split) + '\n'
        cen_lines.append(icoor_cen_line)
        #Centroid version line
        cen_lines.insert(header_index, '# centroid version 0.1\n')

        # Write modified lines to output file
        if name[0]!='7':
            with open(output_centroid_file, 'w') as cen_out:
                cen_out.writelines(cen_lines)


if __name__ == "__main__":
    import sys
    if sys.argv[1] == 'help':
        print("Usage: python peptoid_centroid_gen.py '<input_files>' <output_directory>")
    else:
        fa_input = sys.argv[1]
        centroid_files_dir = sys.argv[2]
        main(fa_input, centroid_files_dir)
        print('Peptoid centroid params file generated, check output file at output directory')
        