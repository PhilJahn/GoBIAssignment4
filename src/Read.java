import java.util.ArrayList;

import AugmentedTree.Interval;
import AugmentedTree.IntervalTree;
import net.sf.samtools.*;

public class Read implements Interval{

	int start;
	int stop;
	String readname;
	//false =forward, true = reverse
	boolean strand;
	
	IntervalTree<RegionBlock> alignment_blocks;
	
	public Read(SAMRecord samRecord1, SAMRecord samRecord2) {
		SAMRecord fop;
		SAMRecord sop;
		if(samRecord1.getFirstOfPairFlag()){
			fop = samRecord1;
			sop = samRecord2;
		}
		else{
			sop = samRecord1;
			fop = samRecord2;
		}
		
		start = Math.min(fop.getAlignmentStart(), sop.getAlignmentStart());
		stop = Math.max(fop.getAlignmentEnd(), sop.getAlignmentEnd());
		readname = fop.getReadName();
		strand = fop.getReadNegativeStrandFlag();
		
		ArrayList<AlignmentBlock> fop_blocks = new ArrayList<AlignmentBlock>(fop.getAlignmentBlocks());
		ArrayList<AlignmentBlock> sop_blocks = new ArrayList<AlignmentBlock>(sop.getAlignmentBlocks());

		
	}

	@Override
	public int getStart() {
		return start;
	}

	@Override
	public int getStop() {
		return stop;
	}

}
