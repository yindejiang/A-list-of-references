
import astropy.io.fits as fits
import sys

def check_nline(fits_filename):
    # Open the FITS file
    with fits.open(fits_filename) as hdulist:
        # Get the header of the first extension
        header1 = hdulist[1].header
        
        # Get the total number of time samples
        nline = header1['NAXIS2']
        
        # Print the valid range
        print(f"Valid time sample index range: 0 to {nline - 1}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python check_nline.py <fits_filename>")
        sys.exit(1)
    
    fits_filename = sys.argv[1]
    check_nline(fits_filename)
