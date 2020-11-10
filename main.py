import sys
import os
import numpy as np
import math
import pandas as pd


def is_diag_code_present(diag_code, list_of_general_codes):
	for general_code in list_of_general_codes:
		if diag_code.startswith(general_code):
			return True
	return False

def get_stage_of_patient(unique_diag_codes, early_diag_codes, late_diag_codes):
	is_at_least_one_code_contained_in_early_diag_codes = False
	is_at_least_one_code_contained_in_late_diag_codes = False
	
	for diag_code in unique_diag_codes:
		# Check if at least one diag_code is contained in early_diag_codes
		if is_diag_code_present(diag_code, early_diag_codes):
			is_at_least_one_code_contained_in_early_diag_codes = True
			break
	
	for diag_code in unique_diag_codes:
		# Check if at least one diag_code is contained in late_diag_codes
		if is_diag_code_present(diag_code, late_diag_codes):
			is_at_least_one_code_contained_in_late_diag_codes = True
			break
	
	if is_at_least_one_code_contained_in_early_diag_codes == True and is_at_least_one_code_contained_in_late_diag_codes == False:
		return 'early'
	
	if is_at_least_one_code_contained_in_early_diag_codes == True and is_at_least_one_code_contained_in_late_diag_codes == True:
		return 'late'
	
	return None

def calc_RR_and_PCC(CO, P, N):
	num_codes = CO.shape[0]
	RR = np.zeros((num_codes,num_codes), dtype=float)
	PCC = np.zeros((num_codes,num_codes), dtype=float)
	
	for i in range(0,num_codes):
		for j in range(0,num_codes):
			RR[i][j] = float(CO[i][j]*N) / float(P[i]*P[j])
			if N-P[i] == 0 or N-P[j] == 0:
				PCC[i][j] = -99999999 # INF
			else:
				PCC[i][j] = float(CO[i][j]*N - P[i]*P[j]) / float( math.sqrt(P[i]) * math.sqrt(P[j]) * math.sqrt(N-P[i]) * math.sqrt(N-P[j]) )
			RR[j][i] = RR[i][j]
			PCC[j][i] = PCC[i][j]
	
	return [RR, PCC]

def get_ranked_list_of_diags(CO, map_index_to_code):
	CO_stage = []
	for j in range(0,CO.shape[1]):
		CO_stage.append( [map_index_to_code[j], CO[0][j]] )
	
	CO_stage = pd.DataFrame(np.array(CO_stage), columns = ['code','freq'])
	CO_stage.freq = CO_stage.freq.astype(int)
	CO_stage = CO_stage.sort_values(by='freq', ascending=False)
	return CO_stage
	
def get_ranked_list_of_procs(P_procs, map_index_to_code):
	CO_stage = []
	for i in range(0,len(P_procs)):
		CO_stage.append( [map_index_to_code[i], P_procs[i]] )

	CO_stage = pd.DataFrame(np.array(CO_stage), columns = ['code','freq'])
	CO_stage.freq = CO_stage.freq.astype(int)
	CO_stage = CO_stage.sort_values(by='freq', ascending=False)
	return CO_stage

def save_mat(mat, map_index_to_code, out_file_path):
	num_codes = mat.shape[0]
	out_file = open(out_file_path, 'w')
	
	for j in range(0,num_codes):
		out_file.write(',' + map_index_to_code[j])
	out_file.write('\n')
	for i in range(0,num_codes):
		out_file.write(map_index_to_code[i])
		for j in range(0,num_codes):
			out_file.write(',' + str(mat[i][j]))
		out_file.write('\n')
	
	out_file.close()

def main():
	dataset_path = sys.argv[1]
	cancer_name = sys.argv[2]
	gender = sys.argv[3]
	stage = sys.argv[4]
	print('dataset_path:', dataset_path)
	print('cancer_name:', cancer_name)
	print('gender:', gender)
	print('stage:', stage)
	
	early_diag_codes = None
	late_diag_codes = None
	if cancer_name == 'colorectal-cancer' and gender in ['male', 'female']:
		early_diag_codes = ['153', '154']
		late_diag_codes = ['196', '197', '198']
	elif cancer_name == 'non-hodkins-lymphoma' and gender == 'male':
		early_diag_codes = ['2028']
		late_diag_codes = ['---invalid---']
	elif cancer_name == 'leukemia' and gender in ['male', 'female']:
		early_diag_codes = ['205', '206', '207', '208']
		late_diag_codes = ['---invalid---']
	elif cancer_name == 'leukemia-limphoid' and gender in ['male', 'female']:
		early_diag_codes = ['204']
		late_diag_codes = ['---invalid---']
	elif cancer_name == 'lung-cancer' and gender in ['male', 'female']:
		early_diag_codes = ['162']
		late_diag_codes = ['196', '197', '198']
	elif cancer_name == 'liver-cancer' and gender in ['male', 'female']:
		early_diag_codes = ['155']
		late_diag_codes = ['196', '197', '198']
	elif cancer_name == 'breast-cancer' and gender == 'female':
		early_diag_codes = ['174']
		late_diag_codes = ['196', '197', '198']
	elif cancer_name == 'thyroid-cancer' and gender in ['male', 'female']:
		early_diag_codes = ['193']
		late_diag_codes = ['196', '197', '198']
	elif cancer_name == 'corpus-uteri' and gender == 'female':
		early_diag_codes = ['182']
		late_diag_codes = ['196', '197', '198']
	elif cancer_name == 'colon-cancer' and gender in ['male', 'female']:
		early_diag_codes = ['153']
		late_diag_codes = ['196', '197', '198']
	elif cancer_name == 'rectal-cancer' and gender in ['male', 'female']:
		early_diag_codes = ['154']
		late_diag_codes = ['196', '197', '198']
	elif cancer_name == 'cervix-uteri-cancer' and gender == 'female':
		early_diag_codes = ['180']
		late_diag_codes = ['196', '197', '198']
	elif cancer_name == 'myeloid-leukemia' and gender in ['male', 'female']:
		early_diag_codes = ['205']	
		late_diag_codes = ['---invalid---']
	elif cancer_name == 'monocytic-leukemia' and gender in ['male', 'female']:
		early_diag_codes = ['206']
		late_diag_codes = ['---invalid---']
	else:
		print('Invalid disease name!')
		print('Exiting ...')
		sys.exit()
	
	input_stage_diag_codes = None
	if stage == 'early':
		input_stage_diag_codes = early_diag_codes
	elif stage == 'late':
		input_stage_diag_codes = late_diag_codes
	else:
		print('Invalid stage!')
		print('Exiting ...')
		sys.exit()
	
	out_dir = 'output/' + cancer_name + '_' + gender + '_' + stage + '/'
	
	if gender == 'male':
		gender = '0'
	elif gender == 'female':
		gender = '1'
	else:
		print('Invalid gender!')
		print('Exiting ...')
		sys.exit()
	
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	
	# ====================================================================================================================
	
	map_patient_id_to_unique_diag_codes = {}
	map_patient_id_to_unique_proc_codes = {}
	
	with open(dataset_path) as file:
		for line in file:
			tokens = line.replace('\n','').split(',')
			patient_id = tokens[50].strip()
			patient_gender = tokens[45].strip()

			if patient_id != '-99999999' and patient_gender == gender:
				
				if patient_id not in map_patient_id_to_unique_diag_codes:
					map_patient_id_to_unique_diag_codes[patient_id] = []
				if patient_id not in map_patient_id_to_unique_proc_codes:
					map_patient_id_to_unique_proc_codes[patient_id] = []
				
				diag_codes = tokens[5:30]
				for code in diag_codes:
					curr_diag_code = code.strip()
					if curr_diag_code.lower() != '\\n':
						if curr_diag_code not in map_patient_id_to_unique_diag_codes[patient_id]:
							map_patient_id_to_unique_diag_codes[patient_id].append(curr_diag_code)
				
				proc_codes = tokens[30:45]
				for code in proc_codes:
					curr_proc_code = code.strip()
					if curr_proc_code.lower() != '\\n':
						if curr_proc_code not in map_patient_id_to_unique_proc_codes[patient_id]:
							map_patient_id_to_unique_proc_codes[patient_id].append(curr_proc_code)
	
	# ====================================================================================================================
	
	# Filter only those patients that have a <stage>-stage <cancer_name>
	map_diag_code_to_index = {}
	map_index_to_diag_code = {}
	map_diag_code_to_index['_'.join(input_stage_diag_codes)] = 0
	map_index_to_diag_code[0] = '_'.join(input_stage_diag_codes)
	diag_index = 1
	
	map_proc_code_to_index = {}
	map_index_to_proc_code = {}
	proc_index = 0
	
	FILTERED_map_patient_id_to_unique_diag_codes = {}
	FILTERED_map_patient_id_to_unique_proc_codes = {}
	
	for patient_id in map_patient_id_to_unique_diag_codes:
		
		curr_patient_stage = get_stage_of_patient(map_patient_id_to_unique_diag_codes[patient_id],
							  early_diag_codes, late_diag_codes)
		
		if curr_patient_stage == stage:
			
			FILTERED_map_patient_id_to_unique_diag_codes[patient_id] = []
			for diag_code in map_patient_id_to_unique_diag_codes[patient_id]:
				if is_diag_code_present(diag_code, input_stage_diag_codes) == False:
					FILTERED_map_patient_id_to_unique_diag_codes[patient_id].append(diag_code)
			
			FILTERED_map_patient_id_to_unique_proc_codes[patient_id] = []
			for proc_code in map_patient_id_to_unique_proc_codes[patient_id]:
				FILTERED_map_patient_id_to_unique_proc_codes[patient_id].append(proc_code)
			
			for diag_code in FILTERED_map_patient_id_to_unique_diag_codes[patient_id]:
				if diag_code not in map_diag_code_to_index:
					map_diag_code_to_index[diag_code] = diag_index
					map_index_to_diag_code[diag_index] = diag_code
					diag_index += 1
	
			for proc_code in FILTERED_map_patient_id_to_unique_proc_codes[patient_id]:
				if proc_code not in map_proc_code_to_index:
					map_proc_code_to_index[proc_code] = proc_index
					map_index_to_proc_code[proc_index] = proc_code
					proc_index += 1		
	
	# ====================================================================================================================
	
	num_diags = len(map_diag_code_to_index)
	CO_diags = np.zeros((num_diags,num_diags), dtype=int)
	P_diags = np.zeros((num_diags,), dtype=int)

	num_procs = len(map_proc_code_to_index)
	CO_procs = np.zeros((num_procs,num_procs), dtype=int)
	P_procs = np.zeros((num_procs,), dtype=int)

	N = len(FILTERED_map_patient_id_to_unique_diag_codes.keys())

	for patient_id in FILTERED_map_patient_id_to_unique_diag_codes:
		
		curr_unique_diag_codes = FILTERED_map_patient_id_to_unique_diag_codes[patient_id]
		P_diags[0] += 1
		for i in range(0,len(curr_unique_diag_codes)):
			index_i = map_diag_code_to_index[curr_unique_diag_codes[i]]
			for j in range(i,len(curr_unique_diag_codes)):
				index_j = map_diag_code_to_index[curr_unique_diag_codes[j]]
				if curr_unique_diag_codes[i] != curr_unique_diag_codes[j]:
					CO_diags[index_i][index_j] += 1
					CO_diags[index_j][index_i] += 1
			CO_diags[0][index_i] += 1
			CO_diags[index_i][0] += 1
			P_diags[index_i] += 1
		
		curr_unique_proc_codes = FILTERED_map_patient_id_to_unique_proc_codes[patient_id]
		for i in range(0,len(curr_unique_proc_codes)):
			index_i = map_proc_code_to_index[curr_unique_proc_codes[i]]
			for j in range(i,len(curr_unique_proc_codes)):
				index_j = map_proc_code_to_index[curr_unique_proc_codes[j]]
				if curr_unique_proc_codes[i] != curr_unique_proc_codes[j]:
					CO_procs[index_i][index_j] += 1
					CO_procs[index_j][index_i] += 1
			P_procs[index_i] += 1
	
	CO_diags_stage = get_ranked_list_of_diags(CO_diags, map_index_to_diag_code)
	CO_procs_stage = get_ranked_list_of_procs(P_procs, map_index_to_proc_code)
	
	# Save N, P_diags, P_procs
	out_file = open(out_dir + 'N.csv', 'w')
	out_file.write(str(N) + '\n')
	out_file.close()
	
	out_file = open(out_dir + 'P_diags.csv', 'w')
	out_file.write('diag_code,prevalence\n')
	for i in range(0,num_diags):
		out_file.write(map_index_to_diag_code[i] + ',' + str(P_diags[i]) + '\n')
	out_file.close()
	
	out_file = open(out_dir + 'P_procs.csv', 'w')
	out_file.write('proc_code,prevalence\n')
	for i in range(0,num_procs):
		out_file.write(map_index_to_proc_code[i] + ',' + str(P_procs[i]) + '\n')
	out_file.close()
	
	# Save the ranked lists of diags/procs for the chosen stage
	out_file = open(out_dir + 'CO_diags_' + stage + '.csv', 'w')
	out_file.write('diag_code,frequency\n')
	for index, row in CO_diags_stage.iterrows():
		out_file.write(str(row['code']) + ',' + str(row['freq']) + '\n')
	out_file.close()
	
	out_file = open(out_dir + 'CO_procs_' + stage + '.csv', 'w')
	out_file.write('proc_code,frequency\n')
	for index, row in CO_procs_stage.iterrows():
		out_file.write(str(row['code']) + ',' + str(row['freq']) + '\n')
	out_file.close()
	
	[RR_diags, PCC_diags] = calc_RR_and_PCC(CO_diags, P_diags, N)
	
	# Save CO_diags, RR_diags, PCC_diags and CO_procs
	save_mat(CO_diags, map_index_to_diag_code, out_dir + 'CO_diags.csv')
	save_mat(RR_diags, map_index_to_diag_code, out_dir + 'RR_diags.csv')
	save_mat(PCC_diags, map_index_to_diag_code, out_dir + 'PCC_diags.csv')
	save_mat(CO_procs, map_index_to_proc_code, out_dir + 'CO_procs.csv')
	print('\nDone!')

if __name__== "__main__":
	main()
