# Automated bot task
## Introduction
ProteinBoxBot is a bot framework which acts as a code base of a family of bots, which aims at making wikidata the central hub of the life sciences. 
To achieve this a family of bots, each with a specific task and domain (i.e. genes, proteins, diseases, and drugs) need to run in frequent update cycle. 
In the genewiki project, a continuous integration tool is used to achieve an orchestrated sequence of automated bot runs, to keep wikidata up-to-date with content
from authoritative resources. This document describe the steps involved to get a specific bot with a specific task to be included inthe automatic bot cycle, maintained by 
the jenkins-instance run at http://genewiki.micel.io:8080.

### Step 1: Create a folder in one of the subfolders of this folder.
This folder groups different bots according to their main topic, being genes, proteins, diseases and drugs. Subsequent grouping can be possible, for example distinguishing between mamalian and microbial genes (i.e. genes/mammals and genes/bacteria)

### Step 2: Adapt your bot code to reflect the above mentioned folder structure.
To allow for testing if new features or improvements on ProteinBoxBot_Core does not intervere with the automated bot cycles run in jenkins, a separate copy of the ProteinBoxBot_Core is kept 
in the automated_bots folder. Bots added to this folder should use this copy of ProteinBoxBot_Core and not the development version.

### Step 3: Enable logs
ProteinBoxBot comes with a log call. The following template can be used to enable logging:

```
try:
  # bot code processing a wikidata item. 
except Exception as e:
	PBB_Core.WDItemEngine.log('ERROR', '{main_data_id}, "{exception_type}", "{message}", {wd_id}, {duration}'.format(
                        main_data_id=diseaseClass.do_id,
                        exception_type=type(e),
                        message=e.__str__(),
                        wd_id='-',
                        duration=time.time() - self.start
                    ))
```

### Step 4: Enable storing of the password of wikidata as an environment variable.
Each time a bot is run a copy of that bot is pulled from bitbucket. This means that passwords are not stored locally as part of config file or other internal settings mechanism, or as a command line variable.  In the current jenkins workflows, passwords are submitted as environment variable. The following python snippet can be used to get the environemnt variable:
```
import os
login = PBB_login.WDLogin("username", os.environ['wikidataApi'])
```

### Step 5: Create a README.md file with contact details and update frequence
The readme file should at least have the following details:
```
# Bot details
** source: ** Name of the original source (e.g. Disease ontology)
** accesspoint: ** Access point where the authoritative content is comming from (e.g. http://purl.obolibrary.org/obo/doid.owl)
** update frequency: ** How often should the bot run (e.g. Every two weeks)
** curator:** Who to contact on content issues?
** bot-operator:** Who to contact on bot issues?
``` 

### Step 6: request launch
Send a note to me (Andra Waagmeester) that a new bot task is ready to run on jenkins in automatic mode. 
