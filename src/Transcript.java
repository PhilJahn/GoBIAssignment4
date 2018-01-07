import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;

import AugmentedTree.IntervalTree;

public class Transcript extends RegionVector{
	
	private IntervalTree<Region> introns; 

	public Transcript(int start, int stop, Annotation annotation)  {
		super(start, stop, annotation);
	}
	
	public int getProtNum(){
		ArrayList<Annotation> annos = new ArrayList<Annotation>();
		for(Region r: this.getRegionsTree()){
			Annotation anno = r.getAnnotation();
			if(!annos.contains(anno)){
				annos.add(anno);
			}
		}
		return annos.size();
		
	}
	
	public boolean setIntrons(){
		if(this.getRegionsTree().size() > 1){
			introns = this.invert().getRegionsTree();
			return true;
		}
		else{
			return false;
		}
	}
	
	public IntervalTree<Region> getIntrons(){
		return introns;
	}
	
	public boolean equals(Object o){
		if(o.getClass() == this.getClass()){
			return equals((Transcript) o);
		}
		return false;
	}
	
	public boolean equals(Transcript t){
		return this.getAnnotation().equals(t.getAnnotation());
	}
	
	public int hashCode(){
		return this.getAnnotation().hashCode();
	}
	
	public String toString(){
		return "Transcript: " + super.toString();
	}

	public boolean inTranscript(Read curRead) {
		
//		System.out.println(curRead.toString());
		
		
		boolean in = true;
		ArrayList<RegionBlock> exonsfop = curRead.getAlignmentBlocksFoP();
		ArrayList<RegionBlock> exonssop = curRead.getAlignmentBlocksSoP();
		ArrayList<RegionBlock> intronsR = curRead.getIntronBlocks();
		IntervalTree<Region> regions = this.getRegionsTree();
		
		if(curRead.getReadName().equals("135")){
			System.out.println("135");
			System.out.println(this.getRegions().toString());
			System.out.println(exonsfop.toString());
			System.out.println(exonssop.toString());
			System.out.println(intronsR.toString());
		}
		
		in = inTranscriptExon(exonsfop);
		if(curRead.getReadName().equals("135")){
			System.out.println("FOP " + in);
		}
		if(in){
			in = inTranscriptExon(exonssop);
			if(curRead.getReadName().equals("135")){
				System.out.println("SOP " + in);
			}
		}
		
		
		if(in && intronsR.size() > 0){
			HashSet<Region> cont;
			for(RegionBlock intron: intronsR){
				cont = regions.getIntervalsSpannedBy(intron.getStart()-1,intron.getStop()-1, new HashSet<Region>());
				if(cont.size() != 0){
					if(curRead.getReadName().equals("135")){
						System.out.println("Intron " + false);
					}
					return false;
				}
			}
			if(curRead.getReadName().equals("135")){
				System.out.println("Intron " + in);
			}
		}
		
		
		return in;
	}
	
	private boolean inTranscriptExon(ArrayList<RegionBlock> exons){
		boolean in = true;
		ArrayList<Region> cont;
		
		IntervalTree<Region> regions = this.getRegionsTree();
		
		RegionBlock curExon = exons.get(0);
		
		if(exons.size() > 1){
			cont = regions.getIntervalsEndAt(curExon.getStop()-1, new ArrayList<Region>());
			if(cont.size() > 0){
				in &= cont.get(0).getStart() <= curExon.getStart();
			}
			else{
				return false;
			}
			for(int i = 1; i < exons.size()-1; i++){
				if(in){
					curExon = exons.get(i);
					cont = regions.getIntervalsEqual(curExon.getStart(),curExon.getStop()-1, new ArrayList<Region>());
					in &= cont.size() > 0;
				}
			}
			if(in){
				curExon = exons.get(exons.size()-1);
				cont = regions.getIntervalsBeginAt(curExon.getStart(), new ArrayList<Region>());
				if(cont.size() > 0){
					in &= cont.get(0).getStop() >= curExon.getStop()-1;
				}
				else{
					return false;
				}
			}
		}
		else{
			cont = regions.getIntervalsSpanning(curExon.getStart(),curExon.getStop()-1, new ArrayList<Region>());
			in &= cont.size() > 0;
		}
		return in;
	}
	
}
