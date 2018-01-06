import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.SortedSet;

import AugmentedTree.Interval;
import AugmentedTree.IntervalTree;
import net.sf.samtools.*;

public class Read implements Interval{

	int start;
	int stop;
	String readname;
	//false =forward, true = reverse
	boolean strand;
	boolean incons;
	
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
		
		this.start = Math.min(fop.getAlignmentStart(), sop.getAlignmentStart());
		int o_start = Math.max(fop.getAlignmentStart(), sop.getAlignmentStart());
		this.stop = Math.max(fop.getAlignmentEnd(), sop.getAlignmentEnd());
		int o_stop = Math.min(fop.getAlignmentEnd(), sop.getAlignmentEnd());
		
		readname = fop.getReadName();
		strand = fop.getReadNegativeStrandFlag();
		
		ArrayList<RegionBlock> exons_fop = new ArrayList<RegionBlock>();
		for (AlignmentBlock block : fop.getAlignmentBlocks()){
			int start = block.getReferenceStart();
			int stop = start + block.getLength();
			exons_fop.add(new RegionBlock(start,stop));
		}
		
//		exons_fop.clear();
//		
//		exons_fop.add(new RegionBlock(114357090, 114347830));
//		exons_fop.add(new RegionBlock(114347890, 114347902));
		
		exons_fop.sort(new RegionBlockComparator());
		
		ArrayList<RegionBlock> introns = new ArrayList<RegionBlock>();
		for(int i =1 ; i < exons_fop.size(); i++){
			int start = exons_fop.get(i-1).getStop();
			int stop = exons_fop.get(i).getStart();
			introns.add(new RegionBlock(start+1,stop));
		}

		ArrayList<RegionBlock> exons_sop = new ArrayList<RegionBlock>();
		for (AlignmentBlock block : sop.getAlignmentBlocks()){
			int start = block.getReferenceStart();
			int stop = start + block.getLength();
			exons_sop.add(new RegionBlock(start,stop));
		}
		
//		exons_sop.clear();
//		
//		exons_sop.add(new RegionBlock(114357090, 114357336));
//		exons_sop.add(new RegionBlock(114357090, 114357168));
//		exons_sop.add(new RegionBlock(114347890, 114347902));
		
		exons_sop.sort(new RegionBlockComparator());
		
		for(int i =1 ; i < exons_sop.size(); i++){
			int start = exons_sop.get(i-1).getStop();
			int stop = exons_sop.get(i).getStart();
			introns.add(new RegionBlock(start+1,stop));
		}
		
		IntervalTree<RegionBlock> exons = new IntervalTree<RegionBlock>();
		exons.addAll(exons_sop);
		exons.addAll(exons_fop);
		
		for( RegionBlock intron : introns){
			if(intron.getStart() <= o_stop && intron.getStop() > o_start){
				
				System.out.println(intron);
				
				HashSet<RegionBlock> cont = new HashSet<RegionBlock>();
				cont = exons.getIntervalsIntersecting(intron.getStart()+1, intron.getStop()-1, cont);
				ArrayList<RegionBlock> bord = new ArrayList<RegionBlock>();
				bord =exons.getIntervalsIntersecting(intron.getStart(), intron.getStop(), bord);
				
				if(cont.size() != 0){
					System.out.println("cont " + cont.size());
					incons = true;
				}
				
				if(bord.size() != 4){
					System.out.println("bord " + bord.size());
					incons = true;
				}
			}
		}
		
		alignment_blocks = new IntervalTree<RegionBlock>();
		
		Iterator<Set<RegionBlock>> exon_it = exons.groupIterator();
		while(exon_it.hasNext()){
			Set<RegionBlock> overlap = exon_it.next();
			RegionBlock[] overlapArray = new RegionBlock[overlap.size()];
			overlapArray = overlap.toArray(overlapArray);
			Arrays.sort(overlapArray,new RegionBlockComparator());
			int start = overlapArray[0].getStart();
			int stop = overlapArray[overlapArray.length-1].getStop();
			
			alignment_blocks.add(new RegionBlock(start,stop));
		}
		if(incons){
			System.out.println("FoP: " + exons_fop.toString());
			System.out.println("SoP: " + exons_sop.toString());
			
			System.out.println("Intron: " + introns.toString());
			
			System.out.println("Ges: "+ this.start + ":" + this.stop);
			System.out.println("Other: " + o_start + ":" + o_stop);
			System.out.println(alignment_blocks.toTreeString());
		}
	}

	@Override
	public int getStart() {
		return start;
	}

	@Override
	public int getStop() {
		return stop;
	}

	class RegionBlockComparator implements Comparator<RegionBlock>
	{
	    public int compare(RegionBlock x1, RegionBlock x2)
	    {
	        return x1.getStart() - x2.getStart();
	    }
	}
	
}
