FROM python:3
MAINTAINER Stefan Verhoeven <s.verhoeven@esciencecenter.nl>

ADD ./ /opt/MARRMoT/

# Install octave and Python deps
RUN apt update && apt install -y octave && pip install -r /opt/MARRMoT/BMI/python/requirements.txt

# Set environment
WORKDIR /data/input
ENV BMI_MODULE=MARRMoTPythonBMI
ENV BMI_CLASS=MARRMoTPythonBMI
ENV BMI_PORT=55555
ENV OCTAVE_MODEL_ROOT=/opt/MARRMoT
ENTRYPOINT ["run-bmi-server", "--path", "/opt/MARRMoT/BMI/python"]
EXPOSE 55555
