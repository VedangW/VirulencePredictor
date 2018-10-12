## About the viruses

The data consists of proteomes from several viruses. The total number of entities for which proteomes are present are 215.

However, some are different versions of the same virus. For instance, take the influenza strain,
```A.HongKong.156.1997.H5N1```.
There are two versions of this virus present in the dataset. These are:
```A.HongKong.156.1997.H5N1.1``` and ```A.HongKong.156.1997.H5N1.2```.

##### Structure for each virus

Each virus is represented by the sequences for its constituent segments, which can be named ```Seg1p1```, ```Seg2p1``` and so on till ```Seg8p2```. These are actually the biological segments, ```PB2```, ```PB1```, ```PB1-F2```, ```PA```, ```HA``` and so on.

These are identified by Identifiers, which can be represented like this:
```>tr|<Type-num>_<Seq-ID>_<Virus-name>_<Seg-num>```

For example, here is a typical identifier from ```A.HongKong.156.1997.H5N1.1```:

```>tr|ASSEM0004_AF046093p1_A/HongKong/156/1997(H5N1)_Seg1p1```

Here the 
```Type-num``` is ```ASSEM0004```,
```Seq-ID``` is ```AF046093p1```,
```Virus-name``` is ```A/HongKong/156/1997/(H5N1)``` and
```Seg-num``` is ```Seg1p1```.

Following are how they vary:
1. For every entity for which proteomes exist, the ```Type-num``` varies in value.
2. For every segment in each entity, the ```Seq-ID``` varies.
It should be noted that the ```Seq-ID``` varies only in the last place if two sequences belong to different types of the same segment, for instance if ```AF046094p1``` is ```Seg2p1``` then ```AF046094p2``` is ```Seg2p2```. Completely different segments would mean that the other parts of the ```Seq-ID``` would also be different. For example, in the above mentioned virus,
```AF046093p1``` is ```Seg1p1``` and ```AF046094p1``` is ```Seg2p1```.
3. Every unique virus has a unique ```Virus-name```. Different versions of the same virus also have the same ```Virus-name```.
4. Each segment in the virus has a different ```Seg-num```. ```Seg-num```s repeat values across all entities.
5. Two identifiers of the same segment for two different versions of the same virus differ in the ```Type-num``` and the ```Seq-ID```.  