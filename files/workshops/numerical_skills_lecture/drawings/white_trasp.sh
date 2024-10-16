for f in *.png ; 
	do ~/magick "$f" -transparent white transp_"$f";
done

