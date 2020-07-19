//get input form the user
dir1 = getDirectory("Choose source Directory ");
list1 = getFileList(dir1);
n= list1.length;
HCcontent=newArray(n);
ECcontent=newArray(n);
label=newArray(n);

for (i=0; i<list1.length; i++) {

		path = dir1+list1[i];
		open(path);
		run("8-bit");
		run("Clear Results");

		setSlice(13);
		//run("Threshold...");
		setAutoThreshold("IsoData dark");
		run("Set Measurements...", "area mean standard stack limit display redirect=None decimal=4");
		setSlice(1);

		run("Analyze Particles...", "size=20-Infinity circularity=0.00-1.00 show=Nothing display exclude stack");
		
		num=nResults;
		nuc_area=newArray(num);
		low_thresh=newArray(num);
		slice_num=newArray(num);
		hcarea=newArray(num);

		for(k=0; k<num; k++){
			nuc_area[k]=getResult("Area",k);
			slice_num[k]=getResult("Slice",k);
			low_thresh[k]=getResult("Mean",k) +1.5*getResult("StdDev",k);
		}
		run("Clear Results");

		for(k=0; k<num; k++){
			setSlice(slice_num[k]);
			//run("Threshold...");
			setThreshold(low_thresh[k], 255);
			run("Measure");
		}

		for(k=0; k<num; k++){
			hcarea[k]=getResult("Area",k);
		}

		label[i]=File.getName(path);
		HCcontent[i]=0;
		ECcontent[i]=0;
		for (j=0; j<nResults; j++){
			HCcontent[i]=HCcontent[i]+hcarea[j];
			ECcontent[i]=ECcontent[i]+(nuc_area[j]-hcarea[j]);
		}
		close();

}

run("Clear Results");
for (k=0; k<list1.length; k++) {
	setResult("Label",k ,label[k]);
	setResult("HCcontent",k, HCcontent[k]);
	setResult("ECcontent",k, ECcontent[k]);
	setResult("HC/EC",k, HCcontent[k]/ECcontent[k]);
	
}
