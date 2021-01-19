#ifndef PearlHeader
#define PearlHeader

class Pearl
{
	public:
		int type;
		int gene_age_;						// Time since this gene was born from a non-coding gene		
		int num_horizontal_transfers_;		// HGT counter (transposon jumps between cells)
		int num_vertical_transfers_;		// VGT counter
		int num_jumps_;						// Transposon jump counter (transposon jumps in a cell)
		
		Pearl();
		virtual Pearl *clone() const=0;
		virtual ~Pearl();
		explicit Pearl(const Pearl &bd);
};

#endif
