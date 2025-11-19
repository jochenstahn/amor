EOS-Service
===========

EOS can be used as histogram service to send images to the Nicos instrument control software. 
For that you need to run it on the amor instrument computer:

```bash
amor-nicos {-vv}
```

The instrument config in Nicos needs to configure a Kafka JustBinItImage instance 
for each histogram that should be used:

```python
hist_yz = device('nicos_sinq.devices.just_bin_it.JustBinItImage',
        description = 'Detector pixel histogram over all times',
        hist_topic = 'AMOR_histograms_YZ',
        data_topic = 'AMOR_detector',
        command_topic = 'AMOR_histCommands',
        brokers = ['linkafka01.psi.ch:9092'],
        unit = 'evts',
        hist_type = '2-D SANSLLB',
        det_width = 446,
        det_height = 64,
        ),
hist_tofz = device('nicos_sinq.devices.just_bin_it.JustBinItImage',
        description = 'Detector time of flight vs. z-pixel histogram over all y-values',
        hist_topic = 'AMOR_histograms_TofZ',
        data_topic = 'AMOR_detector',
        command_topic = 'AMOR_histCommands',
        brokers = ['linkafka01.psi.ch:9092'],
        unit = 'evts',
        hist_type = '2-D SANSLLB',
        det_width = 118,
        det_height = 446,
        ),
```

These images have then to be set in the detector configuration as _images_ items:

```                      
images=['hist_yz', 'hist_tofz'],
```
