### Deployment to shinyapps.io
library(rsconnect)
### Before deploying make sure that your computer is authenticated
### You find the secret and the token on shinyapps.io after logging in
## rsconnect::setAccountInfo(name='sec-explorer',
## token='TOKEN',
## secret='SECRET')
# This path needs to be changed to the local folder containing the app
rsconnect::deployApp('~/projects/HeLa_Shiny_App/')
