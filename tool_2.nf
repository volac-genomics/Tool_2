#!/usr/bin/env nextflow


inputContigs = Channel.fromPath("$baseDir/input_genomes/*.fa", type: 'file')

ch_flatContigs = inputContigs.map { [ it.getBaseName(), it ] }


//////////////////////////////////
//        HELP MESSAGE          //
//////////////////////////////////

def helpMessage() {
    log.info """
        Build a custom database and execute mass screening using ABRicate tool.
        Usage: nextflow run tool_2.nf [options]

        Options:
	--customSeq 'xxx'	
            For example: 12P [default: test]
	--identity 'yy'	
            Minimum DNA %identity [default: 75].
	--coverage 'zz'		
            Minimum DNA %coverage [default: 0].

	Example:
	nextflow run tool_2.nf --customSeq '12P' --identity '90' --coverage '75'
        """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


//////////////////////////////////
//        PARAMETERS            //
//////////////////////////////////

// List of abricate databases
custom_sequences = params.customSeq

// Abricate --minid value
abricate_ID = "${params.identity}"

// Abricate --mincov value
abricate_Cov = "${params.coverage}"


//////////////////////////////////
//   SETUP OUTPUT DIRECTORIES   //
//////////////////////////////////

// Generate formatted dateTime
Date today = new Date()
    YearMonthDay = today.format( "yyyy-MM-dd_HH:mm" )
    OutputDir = "${YearMonthDay}"

// Define output path
RunID_output_path = "$baseDir/output_tool2/${custom_sequences}-i${abricate_ID}-c${abricate_Cov}__${OutputDir}"

// Make output directory
File outputDirectory = new File("${RunID_output_path}")
outputDirectory.mkdirs()


//////////////////////////////////
//    SETUP DATABASE INPUT      //
//////////////////////////////////

Channel.from( custom_sequences ).set{ ch_customSeq }

ch_inputSequences = Channel.fromPath("$baseDir/input_custom_seqs/${custom_sequences}_sequences.fa", type: 'file')


//////////////////////////////////
//     WRITE CONFIG TO FILE     //
//////////////////////////////////

// Write run parameters to file, line by line
configFileName = "${RunID_output_path}/run_parameters.txt"

File parametersFile = new File("${configFileName}")
parametersFile.withWriter{ out ->
  params.each {out.println it }
}

Channel.fromPath( "${configFileName}" ).set{ ch_configFile }

////////////////////////////////////////////////////////////////////////////////////


process version {

        publishDir "${RunID_output_path}/", mode: 'copy'

        output:
        file("abricate_version.txt")
        
        script:                                                                 
        """                                                                     
        abricate --version >> abricate_version.txt
        """ 
}


process make_custom_DB {

        tag "$custom_sequences"

        publishDir "${RunID_output_path}", mode: 'copy', pattern: "*.fa"

	input:
	set sequences, file(fred) from ch_customSeq.combine( ch_inputSequences )
		
	output:
	set sequences, file("*") into ch_customDB
        file(fred)

	script:
	"""
	makeblastdb -in $fred -out $custom_sequences/sequences -dbtype nucl -parse_seqids -hash_index
	"""
}


process massScreening {

        tag "$dataset_id"

        cpus 4

        publishDir "${RunID_output_path}/isolate", mode: 'copy'

        input:
        set dataset_id, file(query), database_name, file(database) from ch_flatContigs.combine( ch_customDB )

        output:
        file "${dataset_id}-${custom_sequences}-i${abricate_ID}-c${abricate_Cov}.tab" into ch_abricateOutput

	script:
        """
        abricate --db $custom_sequences --datadir . --minid '${abricate_ID}' --mincov '${abricate_Cov}' $query > ${dataset_id}-${custom_sequences}-i${abricate_ID}-c${abricate_Cov}.tab
        """
}


process massScreeningSummary {

        cpus 4

        publishDir "${RunID_output_path}/summary", mode: 'copy'

        input:
        file('*') from ch_abricateOutput.toList()

        output:
        file "summary_${custom_sequences}-i${abricate_ID}-c${abricate_Cov}.tab" 

        script:
        """
        abricate --summary * > summary_${custom_sequences}-i${abricate_ID}-c${abricate_Cov}.tab
        """
}


//////////////////////////////////// THE END //////////////////////////////////

