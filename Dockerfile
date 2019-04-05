FROM python:3
LABEL maintainer="Stefan Verhoeven <s.verhoeven@esciencecenter.nl>"

# Install octave 4.4.x, Octave packages from Octave forge
RUN echo deb http://deb.debian.org/debian stretch-backports main > /etc/apt/sources.list.d/stretch-backports.list && \
apt update && apt install -t stretch-backports -y python3-pip octave liboctave-dev libnetcdf-dev && \
octave --eval 'pkg install -verbose -forge netcdf io struct statistics optim' && \
echo 'pkg load optim' >> /usr/share/octave/site/m/startup/octaverc

# Install Python deps
ADD BMI/python/requirements.txt /opt/
RUN pip3 install -r /opt/requirements.txt

ADD ./ /opt/MARRMoT/

# Set environment
WORKDIR /data/input
ENV BMI_MODULE=MARRMoTPythonBMI
ENV BMI_CLASS=MARRMoTPythonBMI
ENV BMI_PORT=55555
ENV OCTAVE_MODEL_ROOT=/opt/MARRMoT
ENTRYPOINT ["run-bmi-server", "--path", "/opt/MARRMoT/BMI/python"]
EXPOSE 55555
