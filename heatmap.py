
__usage__ = """
					python heatmap.py
					--exp <EXPRESSION_FILE>
					--genes <GENES_FILE>
					--out <OUTPUT_FILE>
					--samples <SAMPLE_FILE>
					
					optional:
					--zscore <ACTIVATE_ZSCORE_NORMALIZATION>
					"""


import matplotlib.pyplot as plt
import os, sys, re, math
import numpy as np
import seaborn as sns
from pandas import DataFrame

# --- end of imports --- #

def load_data( data_file ):
	
	data = {}
	with open( data_file, "r" ) as f:
		samples = f.readline().strip().split('\t')[1:]
		vals_per_sample = {}
		for sample in samples:
			vals_per_sample.update( { sample: [] } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			exp_per_gene = {}
			for idx, val in enumerate( parts[1:] ):
				exp_per_gene.update( { samples[ idx ]: float( val ) } )
				vals_per_sample[ samples[ idx ] ].append( float( val ) )
			data.update( { parts[0]: exp_per_gene } )
			line = f.readline()
	return data, samples, vals_per_sample


def calculate_z_scores( values ):
	"""! @brief calculate z-scores per gene """
	
	if sum( values ) > 0:
		avg = np.mean( values )
		std = np.std( values )
		zscores = []
		for val in values:
			zscores.append( ( ( val - avg ) / std ) )
		return zscores
	else:
		return values


def construct_data_output_file( data, candidate_genes, candidate_samples, outputfile, gene_name_mapping_table, sample_order, zscore_state ):
	"""! @brief write expression values of all candidate genes into output file """
	
	datamatrix = []
	genes = []
	with open( outputfile, "w" ) as out:
		new_line = [ "gene" ] + sample_order	#sorted( candidate_samples.keys() )
		tissues = new_line[1:]
		out.write( "\t".join( new_line ) + '\n' )
		for gene in candidate_genes:
			new_line = [ gene ]
			tissue_names = []	#tissue names including sample size
			for tissue in tissues:
				tmp_value = []
				for sample in candidate_samples[ tissue ]:
					try:
						tmp_value.append( data[ gene ][ sample ] )
					except KeyError:
						pass
				if len( tmp_value ) > 0:
					new_line.append( sum( tmp_value ) / len( tmp_value ) )
				else:
					new_line.append( 0 )
				#print (tissue + "(n=" + str( len( tmp_value ) ) + ")")
				tissue_names.append( tissue + "_(n=" + str( len( tmp_value ) ) + ")" )
			if zscore_state:
				zscores = calculate_z_scores( new_line[1:] )
				out.write( "\t".join( [ new_line[0] ] + list(map( str, zscores ) )) + '\n' )
				datamatrix.append( zscores )
			else:
				out.write( "\t".join( [ new_line[0] ] + list(map( str, new_line[1:] ) )) + '\n' )
				datamatrix.append( new_line[1:] )
			try:
				genes.append( gene_name_mapping_table[ gene ] )
			except KeyError:
				genes.append( gene )
	return genes, tissue_names, datamatrix


def construct_heatmap( datamatrix, genes, tissues, heatmap_file, zscore_state ):
	"""! @brief construct heatmap from given data matrix """
	
	print ("number of genes for heatmap construction: " + str( len( genes ) ))
	print ("number of samples for heatmap construction: " + str( len( tissues ) ))
	
	df = DataFrame( datamatrix, index=genes, columns=tissues)
	
	fig, ax = plt.subplots( figsize=( 10, 3 ) )
	
	if zscore_state:
		x = sns.heatmap( 	df, ax=ax, linewidths=0.3, annot=True, fmt='.2f', annot_kws={'fontsize':3}, cbar=True, cmap='seismic',  center=0,	#cool
										cbar_kws= { 'label': 'gene expression z-scores', 'shrink': .5 } 
									)
	else:
		x = sns.heatmap( 	df, ax=ax, linewidths=0.3, annot=True, fmt='.2f', annot_kws={'fontsize':3}, cbar=True, cmap='seismic',  center=0,	#cool
										cbar_kws= { 'label': 'gene expression', 'shrink': .5 } 
									)

	x.figure.axes[-1].yaxis.label.set_size(5)
	x.figure.axes[-1].tick_params(axis='both', which='major', labelsize=5)
	x.figure.axes[-1].tick_params(axis='both', which='minor', labelsize=5)
	
	for idx, gene in enumerate( genes ):
		ax.text( -4, idx+0.6, gene, fontsize=5 )
	
	for idx, tissue in enumerate( tissues ):
		ax.text( idx, len( genes ), tissue, fontsize=5, rotation=90 )	#, rotation=90
	
	ax.set_yticklabels( [], rotation=0, fontsize=2 )
	ax.set_xticklabels( [] , rotation=90, fontsize=3  )
	
	ax.xaxis.set_ticks_position('none')
	ax.yaxis.set_ticks_position('none')
	
	plt.yticks( rotation=0 )
	plt.subplots_adjust( left=0.08, right=0.999, top=0.99, bottom=0.6, wspace=0.1 )
	
	plt.savefig( heatmap_file, dpi=300  )
	plt.savefig( heatmap_file.replace( ".pdf", ".svg" ), dpi=300  )


def load_sample_infos( sample_info_file ):
	"""! @brief load sample info from given file """
	
	candidate_samples = {}
	sample_order = []
	
	with open( sample_info_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			sample_order.append( parts[0] )
			candidate_samples.update( { parts[0]: parts[1].split(',') } )
			line = f.readline()
	
	return candidate_samples, sample_order


def main( arguments ):
	"""! @brief run everything """
	
	data_file = arguments[ arguments.index('--exp')+1 ]
	gene_file = arguments[ arguments.index('--genes')+1 ]
	heatmap_file = arguments[ arguments.index('--out')+1 ]
	sample_info_file = arguments[ arguments.index('--samples')+1 ]
	
	data, samples, vals_per_sample = load_data( data_file )
	
	outputfile = heatmap_file.replace( ".pdf", "plotted_values.txt" )
	
	if '--zscore' in arguments:
		zscore_state = True
	else:
		zscore_state = False
	
	gene_name_mapping_table = {}	#could be used to add gene names
	candidate_genes = []
	with open( gene_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			gene_name_mapping_table.update( { parts[0]: parts[1] } )
			candidate_genes.append( parts[0] )
			line = f.readline()
	
	candidate_samples, sample_order = load_sample_infos( sample_info_file )
	
	genes, tissues, datamatrix = construct_data_output_file( data, candidate_genes, candidate_samples, outputfile, gene_name_mapping_table, sample_order, zscore_state )
	
	construct_heatmap( datamatrix, genes, tissues, heatmap_file, zscore_state )


if '--exp' in sys.argv and '--genes' in sys.argv and '--out' in sys.argv and '--samples' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
