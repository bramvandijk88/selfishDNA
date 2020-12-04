#ifndef PearlHeader
#define PearlHeader

class Pearl
{
	public:
		int type;
		Pearl();
		virtual Pearl *clone() const=0;
		virtual ~Pearl();
		explicit Pearl(const Pearl &bd);
};

#endif
