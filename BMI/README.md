# BMI for Octave

## Dependencies:

### Octave

Install Octave on Ubuntu:

**Option 1:** use the Ubuntu Software Manager, click on the GNU Octave and Install.
However, this may noy be the latest version of Octave (Octave 4.4).

For the latest version use

**Option 2:** From the command line
``` 
sudo apt-get install flatpak
flatpak remote-add --if-not-exists flathub https://flathub.org/repo/flathub.flatpakrepo
flatpak install flathub org.octave.Octave
```
Run:
``` 
flatpak run org.octave.Octave
```

### oct2py

```
pip install oct2py
sudo apt-get install gnuplot
sudo apt-get install gnuplot-x11
```


