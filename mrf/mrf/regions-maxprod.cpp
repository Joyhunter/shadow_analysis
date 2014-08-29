// (C) 2002 Marshall Tappen, MIT AI Lab mtappen@mit.edu

#include <cmath>
#include <stdio.h>
#include "MaxProdBP.h"
#include "assert.h"
int numIterRun;
#include "float.h"
// Some of the GBP code has been disabled here


#ifdef __GNUC__
#define _finite(n) finite(n)
#endif


#define mexPrintf printf
#define mexErrMsgTxt printf
OneNodeCluster::OneNodeCluster()
{
}


int OneNodeCluster::numStates;
FLOATTYPE int_exp(FLOATTYPE val, int exp)
{
  FLOATTYPE tmp=1.0;
  for(int i =0; i < exp; i++)
    tmp *=val;
  return tmp;
}

void getPsiMat(OneNodeCluster & /*cluster*/, FLOATTYPE *&destMatrix, 
               int r, int c, MaxProdBP *mrf, int direction)
{
  int mrfHeight = mrf->getHeight();
  int mrfWidth = mrf->getWidth();
  int numLabels = mrf->getNLabels();

  FLOATTYPE *currMatrix = mrf->getScratchMatrix();

  if(((direction==/*OneNodeCluster::*/UP) &&(r==0)) ||
     ((direction==/*OneNodeCluster::*/DOWN) &&(r==(mrfHeight-1)))||
     ((direction==/*OneNodeCluster::*/LEFT) &&(c==0))||
     ((direction==/*OneNodeCluster::*/RIGHT) &&(c==(mrfWidth-1))))
  {
    for(int i=0; i < numLabels * numLabels; i++)
    {
      currMatrix[i] = 1.0;
    }
  }
  else
  {
    int weight_mod = 1;
    if(mrf->varWeights())
/*    {

      switch(direction)
      {
      case OneNodeCluster::LEFT:
        weight_mod = mrf->getHorizWeight(r,c-1);
        break;
      case OneNodeCluster::RIGHT:
        weight_mod = mrf->getHorizWeight(r,c);
        break;
      case OneNodeCluster::UP:
        weight_mod = mrf->getVertWeight(r-1,c);
        break;
      case OneNodeCluster::DOWN:
        weight_mod = mrf->getVertWeight(r,c);
        break;
      }
    }
        */
        {
                if ( direction == /*OneNodeCluster::*/LEFT ) weight_mod = mrf->getHorizWeight(r,c-1);
                else if ( direction == /*OneNodeCluster::*/RIGHT ) weight_mod = mrf->getHorizWeight(r,c);
                else if ( direction == /*OneNodeCluster::*/UP ) weight_mod = mrf->getVertWeight(r-1,c);
                else    weight_mod = mrf->getVertWeight(r,c);
        
        }
    if(weight_mod!=1)
    {
      for(int i = 0; i < numLabels*numLabels; i++)
      {
        currMatrix[i] = int_exp(mrf->getExpV(i),weight_mod);
      }
      destMatrix = currMatrix;
    }
    else
    {
      destMatrix = mrf->getExpV();
    }
  }
}
  
void initOneNodeMsgMem(OneNodeCluster *nodeArray, FLOATTYPE *memChunk, 
                       const int numNodes, const int msgChunkSize)
{
  FLOATTYPE *currPtr = memChunk;
  OneNodeCluster *currNode = nodeArray;
  for(int i = 0; i < numNodes; i++)
  {

    currNode->receivedMsgs[0] = currPtr; currPtr+=msgChunkSize;
    currNode->receivedMsgs[1] = currPtr; currPtr+=msgChunkSize;
    currNode->receivedMsgs[2] = currPtr; currPtr+=msgChunkSize;
    currNode->receivedMsgs[3] = currPtr; currPtr+=msgChunkSize;
    currNode->nextRoundReceivedMsgs[0] = currPtr; currPtr +=msgChunkSize;
    currNode->nextRoundReceivedMsgs[1] = currPtr; currPtr +=msgChunkSize;
    currNode->nextRoundReceivedMsgs[2] = currPtr; currPtr +=msgChunkSize;
    currNode->nextRoundReceivedMsgs[3] = currPtr; currPtr +=msgChunkSize;
  
    currNode++;
  }
}


void OneNodeCluster::ComputeMsgRight(FLOATTYPE *msgDest, int r, int c, MaxProdBP *mrf)
{
  FLOATTYPE *nodeLeftMsg = receivedMsgs[LEFT],
    *nodeDownMsg = receivedMsgs[DOWN],
    *nodeUpMsg =    receivedMsgs[UP];

  FLOATTYPE *psiMat;

  getPsiMat(*this,psiMat,r,c,mrf,/*OneNodeCluster::*/RIGHT);
  
  FLOATTYPE *cmessage = msgDest,
    total = 0;
  
  for(int rightNodeInd = 0; rightNodeInd < numStates; rightNodeInd++)
  {

    *cmessage = 0;
    for(int leftNodeInd = 0; leftNodeInd < numStates; leftNodeInd++)
    {
      FLOATTYPE tmp =   nodeLeftMsg[leftNodeInd] * 
        nodeUpMsg[leftNodeInd] * 
        nodeDownMsg[leftNodeInd] * 
        localEv[leftNodeInd] *
        psiMat[leftNodeInd * numStates + rightNodeInd];

      if(tmp > *cmessage)
        *cmessage = tmp;

      if ((*cmessage != *cmessage)||!(_finite(*cmessage))||(*cmessage < 1e-10))
    {
      mexPrintf("%f %f %f %f %e\n ",nodeLeftMsg[leftNodeInd],nodeUpMsg[leftNodeInd],nodeDownMsg[leftNodeInd] ,
                psiMat[leftNodeInd * numStates + rightNodeInd], localEv[leftNodeInd]);
      
      mexErrMsgTxt("Break Here\n");

      assert(0);
    } 
    }
//      if (*cmessage < 0.000001)
//        *cmessage=0.000001;
    total += *cmessage;
    cmessage++;
  }

  int errFlag = 0;
  for(int i = 0; i < numStates; i++)
  {
    
    msgDest[i] /= total;
    if (msgDest[i] != msgDest[i])
      errFlag=1;
//      if (msgDest[i] < 0.000001)
//        msgDest[i] = 0.000001;
  }
  if(errFlag)
  {
    mexPrintf("%f |",total);
    for(int i = 0; i < numStates; i++)
    {
      mexPrintf("%f ",msgDest[i]);
    }
    mexErrMsgTxt(" ");
    assert(0);
  }
}

// This means, "Compute the message to send left."

void OneNodeCluster::ComputeMsgLeft(FLOATTYPE *msgDest, int r, int c, MaxProdBP *mrf)
{
  FLOATTYPE *nodeRightMsg = receivedMsgs[RIGHT],
    *nodeDownMsg = receivedMsgs[DOWN],
    *nodeUpMsg =    receivedMsgs[UP];

  FLOATTYPE *psiMat;

  getPsiMat(*this,psiMat,r,c,mrf,/*OneNodeCluster::*/LEFT);
  
  FLOATTYPE *cmessage = msgDest,
    total = 0;
  
  for(int leftNodeInd = 0; leftNodeInd < numStates; leftNodeInd++)
  {

    *cmessage = 0;
    for(int rightNodeInd = 0; rightNodeInd < numStates; rightNodeInd++)
    {
       FLOATTYPE tmp =  nodeRightMsg[rightNodeInd] * 
        nodeUpMsg[rightNodeInd] * 
        nodeDownMsg[rightNodeInd] *
        localEv[rightNodeInd] *
        psiMat[leftNodeInd * numStates + rightNodeInd] ; 
       
       if(tmp > *cmessage)
         *cmessage = tmp;

      if ((*cmessage != *cmessage)||!(_finite(*cmessage)))
      {
        mexPrintf("%f %f %f %f \n",nodeRightMsg[rightNodeInd] ,nodeUpMsg[rightNodeInd],nodeDownMsg[rightNodeInd],
                  psiMat[leftNodeInd * numStates + rightNodeInd]);
        
        mexErrMsgTxt("Break Here\n");
      assert(0);
      }
    }
//        if (*cmessage < 0.000001)
//      *cmessage=0.000001;
    total += *cmessage;
    cmessage++;
  }

  for(int i = 0; i < numStates; i++)
  {
    msgDest[i] /= total;
//      if (msgDest[i] < 0.000001)
//        msgDest[i] = 0.000001;
  }
}

void OneNodeCluster::ComputeMsgUp(FLOATTYPE *msgDest, int r, int c, MaxProdBP *mrf)
{
  FLOATTYPE *nodeRightMsg = receivedMsgs[RIGHT],
    *nodeDownMsg = receivedMsgs[DOWN],
    *nodeLeftMsg =    receivedMsgs[LEFT];
 
  FLOATTYPE *psiMat;

  getPsiMat(*this,psiMat,r,c,mrf,/*OneNodeCluster::*/UP);
  
  FLOATTYPE *cmessage = msgDest,
    total = 0;
  
  for(int upNodeInd = 0; upNodeInd < numStates; upNodeInd++)
  {

    *cmessage = 0;
    for(int downNodeInd = 0; downNodeInd < numStates; downNodeInd++)
    {
      FLOATTYPE tmp = nodeRightMsg[downNodeInd] * 
        nodeLeftMsg[downNodeInd] * 
        nodeDownMsg[downNodeInd] * 
        localEv[downNodeInd] *
        psiMat[upNodeInd * numStates + downNodeInd] ; 

       if(tmp > *cmessage)
         *cmessage = tmp;


      if ((*cmessage != *cmessage)||!_finite(*cmessage))
    {
      mexPrintf("%f %f %f %f\n ",nodeRightMsg[downNodeInd],nodeLeftMsg[downNodeInd],nodeDownMsg[downNodeInd],psiMat[upNodeInd * numStates + downNodeInd]);
      
      mexErrMsgTxt("Break Here\n");
      assert(0);
    } 
    }
//      if (*cmessage < 0.000001)
//        *cmessage=0.000001;
    total += *cmessage;
    cmessage++;
  }

  for(int i = 0; i < numStates; i++)
  {
    msgDest[i] /= total;
//      if (msgDest[i] < 0.000001)
//        msgDest[i] = 0.000001;
  }
}

void OneNodeCluster::ComputeMsgDown(FLOATTYPE *msgDest, int r, int c, MaxProdBP *mrf)
{
  FLOATTYPE *nodeRightMsg = receivedMsgs[RIGHT],
    *nodeUpMsg = receivedMsgs[UP],
    *nodeLeftMsg =    receivedMsgs[LEFT];
 
  FLOATTYPE *psiMat;

  getPsiMat(*this,psiMat,r,c,mrf,/*OneNodeCluster::*/DOWN);
  
  FLOATTYPE *cmessage = msgDest,
    total = 0;
  
  for(int downNodeInd = 0; downNodeInd < numStates; downNodeInd++)
  {

    *cmessage = 0;
    for(int upNodeInd = 0; upNodeInd < numStates; upNodeInd++)
    {
      FLOATTYPE tmp =   nodeRightMsg[upNodeInd] * 
        nodeLeftMsg[upNodeInd] * 
        nodeUpMsg[upNodeInd] * 
        localEv[upNodeInd] *
        psiMat[upNodeInd * numStates + downNodeInd] ; 

       if(tmp > *cmessage)
         *cmessage = tmp;

    }
    if (*cmessage != *cmessage)
    {      printf("Break Here\n");
//      if (*cmessage < 0.000001)
//        *cmessage=0.000001;
      assert(0);
    }
   total += *cmessage;
    cmessage++;
  }

  for(int i = 0; i < numStates; i++)
  {
    msgDest[i] /= total;
//      if (msgDest[i] < 0.000001)
//        msgDest[i] = 0.000001;
  }

}




//  void TwoNodeCluster::deliverMsgs()
//  {
//    float *tmp;
  
//    tmp = nextRoundReceivedMsgs[0];
//    nextRoundReceivedMsgs[0] = receivedMsgs[0];
//    receivedMsgs[0] = tmp;

//    tmp = nextRoundReceivedMsgs[1];
//    nextRoundReceivedMsgs[1] = receivedMsgs[1];
//    receivedMsgs[1] = tmp;


//  }

//  void OneNodeCluster::deliverMsgs()
//  {
//      float *tmp;

//      tmp =nextRoundReceivedMsgs[UP];
//      nextRoundReceivedMsgs[UP] =receivedMsgs[UP];
//      receivedMsgs[UP] = tmp;

//      tmp = nextRoundReceivedMsgs[DOWN];
//      nextRoundReceivedMsgs[DOWN] = receivedMsgs[DOWN];
//      receivedMsgs[DOWN] = tmp;

//      tmp = nextRoundReceivedMsgs[LEFT];
//      nextRoundReceivedMsgs[LEFT] = receivedMsgs[LEFT];
//      receivedMsgs[LEFT] = tmp;

//      tmp = nextRoundReceivedMsgs[RIGHT];
//      nextRoundReceivedMsgs[RIGHT] = receivedMsgs[RIGHT];
//      receivedMsgs[RIGHT] = tmp;


//  }


void OneNodeCluster::deliverMsgs()
{
  const double alpha = 0.8,omalpha = 1-alpha;
  int i;
        for( i = 0; i < numStates; i++)
        {
          receivedMsgs[UP][i] = alpha * receivedMsgs[UP][i] + omalpha * nextRoundReceivedMsgs[UP][i];
        }
        for( i = 0; i < numStates; i++)
        {
          receivedMsgs[DOWN][i] = alpha * receivedMsgs[DOWN][i] + omalpha * nextRoundReceivedMsgs[DOWN][i];
        }
        for( i = 0; i < numStates; i++)
        {
          receivedMsgs[LEFT][i] = alpha * receivedMsgs[LEFT][i] + omalpha * nextRoundReceivedMsgs[LEFT][i];
        }
        for( i = 0; i < numStates; i++)
        {
          receivedMsgs[RIGHT][i] = alpha * receivedMsgs[RIGHT][i] + omalpha * nextRoundReceivedMsgs[RIGHT][i];
        }

}

void OneNodeCluster::getBelief(FLOATTYPE *beliefVec)
{
        FLOATTYPE sum=0;
	int i;
        for(i = 0; i < numStates; i++)
        {
                beliefVec[i] = receivedMsgs[UP][i] * receivedMsgs[DOWN][i] *
                  receivedMsgs[LEFT][i] * receivedMsgs[RIGHT][i] * localEv[i];
                if(!_finite(beliefVec[i]))
                {
                  mexPrintf("SPROBLEM!\n %f %f %f  %f\n",receivedMsgs[UP][i] , receivedMsgs[DOWN][i] ,
                         receivedMsgs[LEFT][i] , receivedMsgs[RIGHT][i] );
                  mexErrMsgTxt(" ");
                  assert(0);

                }
                sum += beliefVec[i];
        }

        for( i = 0; i < numStates; i++)
        {
                beliefVec[i] /= sum;
        }
}


void computeMessagesLeftRight(OneNodeCluster *nodeArray, const int numCols, const int /*numRows*/, const int currRow, const FLOATTYPE alpha, MaxProdBP *mrf)
{
  const int numStates = OneNodeCluster::numStates;
  const FLOATTYPE omalpha = 1.0 - alpha;
  int col, i;
  for(col = 0; col < numCols-1; col++)
  {
    nodeArray[currRow * numCols + col].ComputeMsgRight(nodeArray[currRow * numCols + col+1].nextRoundReceivedMsgs[/*OneNodeCluster::*/LEFT],currRow, col, mrf);
    for(i = 0; i < numStates; i++)
    {
      nodeArray[currRow * numCols + col+1].receivedMsgs[/*OneNodeCluster::*/LEFT][i] = 
        omalpha * nodeArray[currRow * numCols + col+1].receivedMsgs[/*OneNodeCluster::*/LEFT][i] + 
        alpha * nodeArray[currRow * numCols + col+1].nextRoundReceivedMsgs[/*OneNodeCluster::*/LEFT][i];
    }
  } 
  for( col = numCols-1; col > 0; col--)
  {
    nodeArray[currRow * numCols + col].ComputeMsgLeft(nodeArray[currRow * numCols + col-1].nextRoundReceivedMsgs[/*OneNodeCluster::*/RIGHT],currRow, col, mrf);
    for(i = 0; i < numStates; i++)
    {
      nodeArray[currRow * numCols + col-1].receivedMsgs[/*OneNodeCluster::*/RIGHT][i] = 
        omalpha * nodeArray[currRow * numCols + col-1].receivedMsgs[/*OneNodeCluster::*/RIGHT][i] + 
        alpha * nodeArray[currRow * numCols + col-1].nextRoundReceivedMsgs[/*OneNodeCluster::*/RIGHT][i];
    }
  } 

}

void computeMessagesUpDown(OneNodeCluster *nodeArray, const int numCols, const int numRows, const int currCol, const FLOATTYPE alpha, MaxProdBP *mrf)
{
  const int numStates = OneNodeCluster::numStates;
  const FLOATTYPE omalpha = 1.0 - alpha;
  int row;
  for(row = 0; row < numRows-1; row++)
  {
    nodeArray[row * numCols + currCol].ComputeMsgDown(nodeArray[(row+1) * numCols + currCol].nextRoundReceivedMsgs[/*OneNodeCluster::*/UP],row, currCol, mrf);
    for(int i = 0; i < numStates; i++)
    {
      nodeArray[(row+1) * numCols + currCol].receivedMsgs[/*OneNodeCluster::*/UP][i] = 
        omalpha * nodeArray[(row+1) * numCols + currCol].receivedMsgs[/*OneNodeCluster::*/UP][i] + 
        alpha * nodeArray[(row+1) * numCols + currCol].nextRoundReceivedMsgs[/*OneNodeCluster::*/UP][i];
    }
  } 
  for(row = numRows-1; row > 0; row--)
  {
    nodeArray[row * numCols + currCol].ComputeMsgUp(nodeArray[(row-1) * numCols + currCol].nextRoundReceivedMsgs[/*OneNodeCluster::*/DOWN], row, currCol, mrf);
    for(int i = 0; i < numStates; i++)
    {
      nodeArray[(row-1) * numCols + currCol].receivedMsgs[/*OneNodeCluster::*/DOWN][i] = 
        omalpha * nodeArray[(row-1) * numCols + currCol].receivedMsgs[/*OneNodeCluster::*/DOWN][i] + 
        alpha * nodeArray[(row-1) * numCols + currCol].nextRoundReceivedMsgs[/*OneNodeCluster::*/DOWN][i];
    }
  } 

}

