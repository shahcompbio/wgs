def dtypes():
    metrics = {
        'cell_id': 'str',
        'unpaired_mapped_reads':'int64',
        'paired_mapped_reads':'int64',
        'unpaired_duplicate_reads':'int64',
        'paired_duplicate_reads':'int64',
        'unmapped_reads':'int64',
        'percent_duplicate_reads':'float64',
        'estimated_library_size':'int64',
        'total_reads':'int64',
        'total_mapped_reads':'int64',
        'total_duplicate_reads':'int64',
        'total_properly_paired':'int64',
        'coverage_breadth':'float64',
        'coverage_depth':'float64',
        
    }

    insert_metrics = {
        'median_insert_size':'int64',
        'mean_insert_size':'float64',
        'standard_deviation_insert_size':'float64'
    }


    dtypes = locals()

    return dtypes
