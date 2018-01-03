import java.awt.List;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

import net.sf.samtools.*;

public class Psi {

	public static void main(String[] args) {
//		long startTime = System.currentTimeMillis();		
		String gtfPath ="";
		String bamPath ="";
		String outputPath ="";
		for(int i =0; i < args.length-1; i++){
			if(args[i].equals("-gtf")){
				gtfPath = args[i+1];
				i++;
			}
			else if(args[i].equals("-bam")){
				bamPath = args[i+1];
				i++; 
			}
			else if(args[i].equals("-o")){
				outputPath = args[i+1];
				i++; 
			}
		}
		
		if(gtfPath.equals("") || bamPath.equals("") || outputPath.equals("")){
			System.out.println("Usage Info:\n-gtf <filepath for GTF>\n-bam <filepath for BAM>\n-o <filepath for ouput>");
		}
		else{
			
			String baiPath = bamPath.concat(".bai");
			
			Path gtfFilePath = Paths.get(gtfPath);
			Path bamFilePath = Paths.get(bamPath);
			Path baiFilePath = Paths.get(baiPath);
			Path outputFilePath = Paths.get(outputPath);
			
			File gtfFile = gtfFilePath.toFile();
			File bamFile = bamFilePath.toFile();
			File baiFile = baiFilePath.toFile();
			
			Psi psi = new Psi(gtfFile, bamFile, baiFile);
		
			String output = psi.getOutput();
			ArrayList<String> outputAsList = new ArrayList<String>();
			outputAsList.add(output);
			try {
				Files.write(outputFilePath, outputAsList, Charset.forName("UTF-8"));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
//		long stopTime = System.currentTimeMillis();
//		System.out.println("Input:" + (stopTime-startTime));	
	}

	private HashMap<Integer,Gene> geneSet;
	
	private SAMFileReader bam_reader;
	
	public Psi (File gtfFile, File bamFile, File baiFile) {
//		long startTime = System.currentTimeMillis();
		
		bam_reader = new SAMFileReader(bamFile,baiFile);
		bam_reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
		
		geneSet = new HashMap<Integer,Gene>();
	    try {
	        BufferedReader br = new BufferedReader (new FileReader(gtfFile));
	        String line;
	        Gene curGene = null;
	        Annotation curGAnno = null;
	        Transcript curTrans = null;
	        Annotation curTAnno = null;
	        Annotation curCAnno = null;
	        ArrayList<Transcript> unknownGene = new ArrayList<Transcript>();
	        ArrayList<Region> unknownTranscript = new ArrayList<Region>();
	        HashMap<String,Annotation> mapAnno = new HashMap<String,Annotation>(); 
	        while ((line = br.readLine()) != null){
	        	String[] lineSplit = line.split("\t");
	        	if(lineSplit.length >= 8){
	        	String[] attrSplit = lineSplit[8].split(";");
	        	HashMap<String,String> attr;

	        	if(lineSplit[2].equals("CDS")){
	        		attr = getAttributes(attrSplit);
	        		String super_super_id = attr.get("gene_id");
	        		String gene_name = attr.get("gene_name");
	        		String super_id = attr.get("transcript_id");
	        		String id = attr.get("protein_id");
	        		String type = lineSplit[1];
	        		char strand = lineSplit[6].charAt(0);
	        		String chr = lineSplit[0];
	        		int start = Integer.parseInt(lineSplit[3]);
	        		int stop = Integer.parseInt(lineSplit[4]);
	        		Region cds;
	        		Annotation cdsAnno;
	        		if(curCAnno != null && curCAnno.getId().equals(id) && curCAnno.getSuperId().equals(super_id) && curCAnno.getSuperSuperId().equals(super_super_id)){
	        			cdsAnno = curCAnno;
	        		}
	        		else{
	        			cdsAnno = mapAnno.get(id+super_id+super_super_id);
	        			if(!(cdsAnno != null)){
	        				cdsAnno = new Annotation(id,chr,strand,super_id,super_super_id,type,gene_name);
	        				mapAnno.put(id+super_id+super_super_id,cdsAnno);
	        			}
	        		}
	        		cds = new Region(start,stop,cdsAnno);
	        		
	        		
	        		if(cdsAnno.isSub(curTAnno)){
	        			curTrans.add(cds);
	        		}
	        		else {
		        		String transcript_name = attr.get("transcript_name");	        			
		        		Annotation transAnno = new Annotation(super_id,transcript_name,chr,strand,super_super_id,type,gene_name);
		        		Transcript trans = new Transcript(start,stop,transAnno);

		        		if(transAnno.isSub(curGAnno)){
		        			curGene.add(trans);
		        		}
		        		else{
			        		Annotation geneAnno = new Annotation(super_super_id,gene_name,chr,strand,type);
			        		Gene gene = new Gene(start, stop, geneAnno);
			        		geneSet.put(gene.hashCode(),gene);
			        		curGene = gene;
			        		curGAnno = geneAnno;
			        		curGene.add(trans);
		        		}
		        		
		        		curTAnno = transAnno;
		        		curTrans = trans;
		        		curTrans.add(cds);
	        		}
	        		
	        	}
	        	else{}
	        	}
	        }
	        br.close();	
	    	
	    } catch (Exception e) {
			e.printStackTrace();
		}
//		long stopTime = System.currentTimeMillis();
//		System.out.println("Input:" + (stopTime-startTime));
	}
	
	public static HashMap<String,String> getAttributes (String[] attrs){
		HashMap <String,String> attrMap = new HashMap <String,String>();
		for(int i = 0; i < attrs.length; i++){
			String curattr = attrs[i];
			curattr = curattr.trim();
			String[] attr = curattr.split("\\s+");
			attr[0] = attr[0].trim();
			attr[1] = attr[1].trim();
			attr[1] = attr[1].replaceAll("\"", "");
			attrMap.put(attr[0],attr[1]);
		}
		return attrMap;
	}
	
	
	public String getOutput(){
//		long startTime = System.currentTimeMillis();
		String tab = "\t";
		String brk = "\n";
		char sep = ':';
		char sip = '|';
		char lin = '-';
//		StringBuilder resultBuilder = new StringBuilder("id\tsymbol\tchr\tstrand\tnprots\tntrans\tSV\tWT\tWT_exons\tWT_prots\tSV_prots\tmin_skipped_exon\tmax_skipped_exon\tmin_skipped_bases\tmax_skipped_bases\n");

		StringBuilder resultBuilder = new StringBuilder("gene\texon\tnum_incl_reads\tnum_excl_reads\tnum_total_reads\tpsi\n");

		String keyString = "ENSG00000158109.10";
		int keyInt = keyString.hashCode();
		Gene gene = geneSet.get(keyInt);
		geneSet.clear();
		geneSet.put(keyInt, gene);
		
		for (Integer key: geneSet.keySet()) {
			Gene curGene = geneSet.get(key);
			Annotation curGAnno = curGene.getAnnotation();
//			System.out.println("getExonSkips");
			ArrayList<ExonSkip> skips = curGene.getExonSkips();
			
//			printExonSkip(resultBuilder, curGene, curGAnno, skips, tab, brk, sep, sip);
			
			String geneid = curGAnno.getId();
			String chr = curGAnno.getChromosome();
			
			for(ExonSkip skip: skips){
				ArrayList<Region> exons = new ArrayList<Region>(skip.getWTExons());
				exons.sort(new StartRegionComparator());
				for( Region exon : exons){
					resultBuilder.append(geneid);
					resultBuilder.append(tab);
					resultBuilder.append(exon.getStart());
					resultBuilder.append(lin);
					resultBuilder.append(exon.getStop()+1);
					
					Iterator<SAMRecord> exon_iterator = bam_reader.queryOverlapping(chr, exon.getStart(), exon.getStop());
					
					int c = 0;
					int d = 0;
					
					HashSet<String> ids = new HashSet<String>();
					
					while(exon_iterator.hasNext()){
						SAMRecord curRec = exon_iterator.next();
						if(ids.add(curRec.getReadName())){
							c++;
							System.out.print(curRec.getReadName() + "\t");
							java.util.List<AlignmentBlock> blocks = curRec.getAlignmentBlocks();
							Iterator<AlignmentBlock> blockIterator = blocks.iterator();
							int blockStart = 0;
							if(blockIterator.hasNext()){
								blockStart = blockIterator.next().getReferenceStart();
								System.out.print(blockStart);
								if(blockStart >= exon.getStart() && blockStart <= exon.getStop()){
									d++;
								}
							}
							while(blockIterator.hasNext()){
								blockStart = blockIterator.next().getReferenceStart();
								System.out.print("|" + blockStart);
								if(blockStart >= exon.getStart() && blockStart <= exon.getStop()){
									d++;
								}
							}
							System.out.println("");
						}
						else{
							
						}
					}
					
					
					
					resultBuilder.append(tab);
					resultBuilder.append(d);
					
					int e = c-d;
					resultBuilder.append(tab);
					resultBuilder.append(e);
					
					resultBuilder.append(tab);
					resultBuilder.append(c);
					
//					int f = 0;
//					Iterator<SAMRecord>  exon_iterator_in = bam_reader.query(chr, exon.getStart(), exon.getStop()+1, true);
//					
//					ids.clear();
//					
//					while(exon_iterator_in.hasNext()){
//						SAMRecord curRec = exon_iterator_in.next();
//						if(ids.add(curRec.getReadName())){
//							f++;
//						}
//					}
					
//					resultBuilder.append(tab);
//					resultBuilder.append(f);
					
					resultBuilder.append(brk);
				}
				
			}
			
		}
//		long stopTime = System.currentTimeMillis();
//		System.out.println("Skips:" + (stopTime-startTime));
		return resultBuilder.toString();
	}
	
	public void printExonSkip(StringBuilder resultBuilder, Gene curGene, Annotation curGAnno, ArrayList<ExonSkip> skips, String tab, String brk, char sep, char sip){

		//		System.out.println("gotExonSkips");
		int tn = curGene.getTranscriptNumber();
		int pn = curGene.getProteinNumber();
		
		String geneid = curGAnno.getId();
		String genename = curGAnno.getName();
		String chr = curGAnno.getChromosome();
		char str = curGAnno.getStrand();
		
		StringBuilder geneInfo = new StringBuilder(geneid);
		geneInfo.append(tab);
		geneInfo.append(genename);
		geneInfo.append(tab);
		geneInfo.append(chr);
		geneInfo.append(tab);
		geneInfo.append(str);
		geneInfo.append(tab);
		geneInfo.append(tn);
		geneInfo.append(tab);
		geneInfo.append(pn);
		geneInfo.append(tab);
		
		for(ExonSkip skip: skips){
			resultBuilder.append(geneInfo);
			resultBuilder.append(skip.getStart());
			resultBuilder.append(sep);
			resultBuilder.append(skip.getStop());
			resultBuilder.append(tab);
			
			ArrayList<Region> introns = new ArrayList<Region>(skip.getWTIntrons());
			introns.sort(new StartRegionComparator());
			resultBuilder.append(introns.get(0).getStart());
			resultBuilder.append(sep);
			resultBuilder.append(introns.get(0).getStop());
			for(int i = 1; i < introns.size(); i++ ){
				Region curIntron = introns.get(i);
				resultBuilder.append(sip);
				resultBuilder.append(curIntron.getStart());
				resultBuilder.append(sep);
				resultBuilder.append(curIntron.getStop());
			}
			resultBuilder.append(tab);
			
			ArrayList<Region> exons = new ArrayList<Region>(skip.getWTExons());
			exons.sort(new StartRegionComparator());
			resultBuilder.append(exons.get(0).getStart());
			resultBuilder.append(sep);
			resultBuilder.append(exons.get(0).getStop());
			for(int i = 1; i < exons.size(); i++ ){
				Region curExon = exons.get(i);
				resultBuilder.append(sip);
				resultBuilder.append(curExon.getStart());
				resultBuilder.append(sep);
				resultBuilder.append(curExon.getStop());
			}
			resultBuilder.append(tab);
			
			ArrayList<String> wt = new ArrayList<String>(skip.getWTProt());
			resultBuilder.append(wt.get(0));
			for(int i = 1; i < wt.size(); i++ ){
				resultBuilder.append(sip);					
				resultBuilder.append(wt.get(i));
			}
			resultBuilder.append(tab);
			
			ArrayList<String> sv = new ArrayList<String>(skip.getSVProt());
			resultBuilder.append(sv.get(0));
			for(int i = 1; i < sv.size(); i++ ){
				resultBuilder.append(sip);					
				resultBuilder.append(sv.get(i));
			}
			resultBuilder.append(tab);
			
			resultBuilder.append(skip.getMinEx());
			resultBuilder.append(tab);
			resultBuilder.append(skip.getMaxEx());
			resultBuilder.append(tab);
			resultBuilder.append(skip.getMinBase());
			resultBuilder.append(tab);
			resultBuilder.append(skip.getMaxBase());
			resultBuilder.append(brk);
		
		}
	}
	
	public ArrayList<Gene> getGenes(){
		ArrayList<Gene> result = new ArrayList<Gene>();
		for (Integer key: geneSet.keySet()) {
		   result.add(geneSet.get(key));
		}
		return result;
	}
	
	class StartRegionComparator implements Comparator<Region>
	{
	    public int compare(Region x1, Region x2)
	    {
	        return x1.getStart() - x2.getStart();
	    }
	}
	

}
