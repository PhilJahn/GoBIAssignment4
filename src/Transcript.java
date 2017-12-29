import java.util.ArrayList;
import java.util.Collection;

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
	
}
