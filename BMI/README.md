# BMI for Octave

This directory contains the [BMI interface](https://bmi-spec.readthedocs.io) of the model.

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


# Building grpc4bmi Docker image

Building the Docker image with [grpc4bmi wrapper](https://grpc4bmi.readthedocs.io) for MARRMoT.

```bash
docker build -t sverhoeven/marrmot .
```

See https://grpc4bmi.readthedocs.io/en/latest/container/usage.html how to use the Docker image.

