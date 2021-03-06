# Visualization-R
     ####Visualization_project.R
     ####irinamahmudjanova
     #### Mon Jan 16 21:41:05 2017
#### The interest of the script is show how to recreate plot from "The Economist" with correlation between corruption and development http://www.economist.com/node/21541178 using the ggplot2 in R.

![economist1](https://cloud.githubusercontent.com/assets/16123495/22446702/db68941e-e702-11e6-861b-8624e58a8483.png)
 
     library(ggplot2)
     library(ggthemes)
     library(data.table)

    ## Import the ggplot2 data.table libraries and use fread to load the csv file 
    ## 'Economist_Assignment_Data.csv' into a dataframe called df (Hint: use drop=1 to skip the first column)
    dc <- read.csv("Economist_Data.csv", header = TRUE, sep = ',', stringsAsFactors = TRUE)  
    dc[ ,1  ] <- NULL
    ## or we can do second way where we will read the data without the first colomn
    dc <- fread("Economist_Data.csv", drop=1)
    ## Check the head of dataframe dc
    ##dc %>% head

    ## Use ggplot() + geom_point() to create a scatter plot object called pl. 
    ## You will need to specify x=CPI and y=HDI and color=Region as aesthetics
    pl <- ggplot(dc, aes(x=CPI, y=HDI, color=Region)) + geom_point()
    print(pl)
![pl1_vis](https://cloud.githubusercontent.com/assets/16123495/22487231/6e04ac72-e7c1-11e6-95b9-7e1686333166.png)
 
#### the same plot we can get if we use color in geom_point
    pl <- ggplot(dc, aes(x=CPI, y=HDI)) + geom_point(aes( color=factor(Region)))
    pl
   
#### we can perform as shown down the same plot changing the size of CPI but we dont need it right now
    pl1 <- pl + geom_point(aes(x=CPI,  y=HDI, color= factor(Region), size = CPI))
    pl1
![pl1_1_vis](https://cloud.githubusercontent.com/assets/16123495/22528352/c6e7df3c-e887-11e6-924c-9f7a3d06889a.png)
 
 We should change the points to be larger empty circles. Let's go back and add arguments to geom_point() and reassign it to pl and try to to figure out what shape and size

     pl2<-ggplot(dc, aes(x=CPI, y=HDI, color=Region)) + geom_point(size =3.5, shape=1)
     pl2
![pl2_vis](https://cloud.githubusercontent.com/assets/16123495/22528358/c947e3c6-e887-11e6-9a82-85fce9df6c04.png)

     pl3 <- pl2+geom_smooth(aes(group=1))
     pl3
![pl3_vis](https://cloud.githubusercontent.com/assets/16123495/22528542/a4cd8aea-e888-11e6-98fa-59e4a74d4b46.png)
we can remove use se=F the gray area and using method linear reg make log line
     
     pl4 <-pl2 + geom_smooth(aes(group=1), method = 'lm', formula=y~log(x), se=F, color='red')
     pl4
![pl4_vis](https://cloud.githubusercontent.com/assets/16123495/22528360/cd3f3bf0-e887-11e6-9918-5d28f00dd456.png)

It's really starting to look similar! But we still need to add labels, we can 
use geom_text! Add geom_text(aes(label=Country)) to pl4 and see what happens. 
Let's take shrink labels, there are way too many of them.

    pl6 <- pl4 + geom_text(aes(label= Region))
    pl6
![pl6_vis](https://cloud.githubusercontent.com/assets/16123495/22528365/d18bbd5a-e887-11e6-8b0a-9ad9ed9388b2.png)

Labeling a subset is actually pretty tricky! So we're just going to give you the
answer since it would require manually selecting the subset of countries we 
want to label!

      pointsToLabel <- c("Russia", "Venezuela", "Iraq", "Myanmar", "Sudan",
                      "Afghanistan", "Congo", "Greece", "Argentina", "Brazil","Morocco",
                      "India", "Italy", "China", "South Africa", "Spane", "Poland",
                      "Botswana", "Kazakhstan", "Rwanda", "France","Switzerland","Romania",
                      "United States", "Germany", "Britain", "Barbados", "Norway", "Japan",
                      "New Zealand", "Singapore", "Uzbekistan", "Kyrgyzstan","Turkey")
     pl7 <- pl4 + geom_text(aes(label = Country), color = "gray20", 
                       data = subset(dc, Country %in% pointsToLabel),check_overlap = TRUE)
    pl7
    
![pl7_vis](https://cloud.githubusercontent.com/assets/16123495/22528372/db6c5924-e887-11e6-9714-993d20ec69cf.png)

Almost there! Still not perfect, but good enough :)
Now let's just add some labels and a theme, set the x and y scales and we're done!
Add theme_bw() to your plot and save this to pl7
    
    pl8 <- pl7+theme_bw()
    pl8


    pl9 <- pl8 + scale_x_continuous(name = "Corruption Perception Index, 2011, (10=least corrupt)", breaks = 1:10)
    pl9


    pl9 <- pl9 + scale_y_continuous(name = "Human Development Business, 2011 (1=Best)", breaks = 1:10)
    pl9
![pl9_vis](https://cloud.githubusercontent.com/assets/16123495/22528391/ec4d680a-e887-11e6-8027-894beaf0f95b.png)
Finally use ggtitle() to add a string as a title.
 
    pl10 <- pl9 + ggtitle("Corruption and Human Development")
    pl10
 ![pl10_vis](https://cloud.githubusercontent.com/assets/16123495/22528398/ef933986-e887-11e6-9a11-2832dc6ea703.png)

Of coerse this plot now absolute copy of posted in The Economist, but it's almost copy and you always can find out how improve any plot in R using wonderful packege ggplot2 :) 
     
     #end
